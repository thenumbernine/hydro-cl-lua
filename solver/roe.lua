local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local FiniteVolumeSolver = require 'solver.fvsolver'

local xNames = table{'x', 'y', 'z'}

-- this can be put in app.lua
local errorType = 'error_t'
local errorTypeCode = [[
typedef struct { 
	real ortho, flux;
} ]]..errorType..';'

local Roe = class(FiniteVolumeSolver)
Roe.name = 'Roe'

-- enable these to verify accuracy
-- disable these to save on allocation / speed
Roe.checkFluxError = false
Roe.checkOrthoError = false

--[[
args specific to Roe:
	checkFluxError
	checkOrthoError
--]]
function Roe:init(args)
	Roe.super.init(self, args)

	if args.checkFluxError ~= nil then
		self.checkFluxError = args.checkFluxError
	end
	if args.checkOrthoError ~= nil then
		self.checkOrthoError = args.checkOrthoError
	end
end

function Roe:createBuffers()
	Roe.super.createBuffers(self)

	-- to get sizeof
	ffi.cdef(self.eqn:getEigenTypeCode())
	ffi.cdef(errorTypeCode)
	
	local realSize = ffi.sizeof(self.app.real)

	self:clalloc('eigenBuf', self.volume * self.dim * ffi.sizeof(self.eqn.eigen_t))
	
	-- debug only
	if self.checkFluxError or self.checkOrthoError then
		local errorTypeSize = ffi.sizeof(errorType)
		self:clalloc('errorBuf', self.volume * self.dim * errorTypeSize)
	end
end

function Roe:createCodePrefix()
	Roe.super.createCodePrefix(self)
	
	self.codePrefix = table{	
		self.codePrefix,
		errorTypeCode,
	}:concat'\n'
end

function Roe:getSolverCode()
	return table{
		Roe.super.getSolverCode(self),
	
		-- before this went above solver/plm.cl, now it's going after it ...
		template(file['solver/roe.cl'], {
			solver = self,
			eqn = self.eqn,
			clnumber = require 'cl.obj.number',
		}),
	}:concat'\n'
end

-- all these are found eqn's cl code
function Roe:refreshSolverProgram()
	Roe.super.refreshSolverProgram(self)

	self.calcEigenBasisKernelObj = self.solverProgramObj:kernel(
		'calcEigenBasis',
		self.eigenBuf,
		self.getULRBuf)

	self.calcFluxKernelObj = self.solverProgramObj:kernel(
		'calcFlux',
		self.fluxBuf,
		self.getULRBuf,
		self.eigenBuf)

	-- TODO put this in solver/solver.lua ?
	if self.eqn.useSourceTerm then
		self.addSourceKernelObj = self.solverProgramObj:kernel{name='addSource', domain=self.domainWithoutBorder}
		self.addSourceKernelObj.obj:setArg(1, self.UBuf)
	end

	if self.checkFluxError or self.checkOrthoError then
		self.calcErrorsKernelObj = self.solverProgramObj:kernel(
			'calcErrors',
			self.errorBuf,
			self.eigenBuf)
	end
end

function Roe:addDisplayVars()
	Roe.super.addDisplayVars(self)

--[=[ TODO add eigen calc code here
	for j,xj in ipairs(xNames) do
		self:addDisplayVarGroup{
			name = 'wave '..xj,
			bufferField = 'waveBuf',
			codePrefix = [[
	int indexInt = ]]..(j-1)..[[ + dim * index;
	const global real* wave = buf + indexInt * numWaves;
]],
			vars = range(0, self.eqn.numWaves-1):map(function(i)
				return {[''..i] = '*value = wave['..i..'];'}
			end),
		}
	end
--]=]

	-- TODO rename to 'getEigenDisplayVarDescs()'
	local eigenDisplayVars = self.eqn:getEigenDisplayVars()
	if eigenDisplayVars and #eigenDisplayVars > 0 then
		for j,xj in ipairs(xNames) do
			self:addDisplayVarGroup{
				name = 'eigen '..xj,
				bufferField = 'eigenBuf',
				type = self.eqn.eigen_t,
				codePrefix = [[
	int indexInt = ]]..(j-1)..[[ + dim * index;
	const global ]]..self.eqn.eigen_t..[[* eigen = buf + indexInt;
]],
				vars = table.map(eigenDisplayVars, function(kv)
					return table.map(kv, function(v,k)
						if k == 'type' then return v, k end
						return v, k
					end)
				end),
			}
		end
	end

	-- TODO add kernels for each side
	if self.checkFluxError or self.checkOrthoError then	
		for j,xj in ipairs(xNames) do
			self:addDisplayVarGroup{
				name = 'error '..xj,
				bufferField = 'errorBuf',
				codePrefix = [[
	int indexInt = ]]..(j-1)..[[ + dim * index;
]],
				useLog = true,
				type = 'error_t',
				vars = {
					{ortho = '*value = buf[indexInt].ortho;'},
					{flux = '*value = buf[indexInt].flux;'},
				},
			}
		end
	end
end

local realptr = ffi.new'real[1]'
local function real(x)
	realptr[0] = x
	return realptr
end

-- NOTICE this adds the contents of derivBuf and does not clear it
function Roe:calcDeriv(derivBuf, dt)
	local dtArg = real(dt)
	
	self:boundary()
	
	if self.usePLM then
		self.calcLRKernelObj.obj:setArg(2, dtArg)
		self.calcLRKernelObj()
	end


	self.calcEigenBasisKernelObj()

	if self.checkFluxError or self.checkOrthoError then
		self.calcErrorsKernelObj()
	end

	self.calcFluxKernelObj.obj:setArg(3, dtArg)
	self.calcFluxKernelObj()

-- [=[ this is from the 2017 Zingale book
	if self.useCTU then
		-- if we're using CTU then ...
		-- 1) calc fluxes based on a slope-limiter method (PLM, etc)
		-- 2) at each interface, integrate each dimension's LR states by all other dimensions' fluxes with a timestep of -dt/2
		--	( don't use the deriv buf because it already has the sum of all dimensions' flux differences)
		self.updateCTUKernelObj.obj:setArg(2, dtArg)
		self.updateCTUKernelObj()

		-- now we need to calcBounds on the ULR
		-- TODO this will break for mirror conditions
		-- because I haven't got the boundary code flexible enough to operate on specific fields within the L & R fields of the ULRBuf
		self.lrBoundaryKernelObj()

		-- 3) use the final LR states to calculate the flux ...

		-- the rest of this matches above
		-- maybe use 'repeat'?
		
		self.calcEigenBasisKernelObj()
		self.calcFluxKernelObj()
	end
--]=]
	
	self.calcDerivFromFluxKernelObj(derivBuf)

	-- addSource adds to the derivative buffer
	if self.eqn.useSourceTerm then
		self.addSourceKernelObj(derivBuf)
	end
end

return Roe

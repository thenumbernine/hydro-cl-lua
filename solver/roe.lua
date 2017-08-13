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

function Roe:createBuffers()
	Roe.super.createBuffers(self)

	-- to get sizeof
	ffi.cdef(self.eqn:getEigenTypeCode())
	ffi.cdef(errorTypeCode)
	
	local realSize = ffi.sizeof(self.app.real)

	self:clalloc('waveBuf', self.volume * self.dim * self.eqn.numWaves * realSize)
	self:clalloc('eigenBuf', self.volume * self.dim * ffi.sizeof(self.eqn.eigen_t))
	self:clalloc('deltaUEigBuf', self.volume * self.dim * self.eqn.numWaves * realSize)
	if self.fluxLimiter[0] > 0 then
		self:clalloc('rEigBuf', self.volume * self.dim * self.eqn.numWaves * realSize)
	end
	
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
		template(file['solver/roe.cl'], {solver=self, eqn=self.eqn}),
	}:concat'\n'
end

-- all these are found eqn's cl code
function Roe:refreshSolverProgram()
	Roe.super.refreshSolverProgram(self)

	self.calcEigenBasisKernel = self.solverProgram:kernel(
		'calcEigenBasis',
		self.waveBuf,
		self.eigenBuf,
		self.getULRBuf)
		
	self.calcDeltaUEigKernel = self.solverProgram:kernel(
		'calcDeltaUEig',
		self.deltaUEigBuf,
		self.getULRBuf,
		self.eigenBuf)

	if self.fluxLimiter[0] > 0 then
		self.calcREigKernel = self.solverProgram:kernel(
			'calcREig',
			self.rEigBuf,
			self.deltaUEigBuf,
			self.waveBuf)
	end
	
	self.calcFluxKernel = self.solverProgram:kernel(
		'calcFlux',
		self.fluxBuf,
		self.getULRBuf,
		self.waveBuf,
		self.eigenBuf,
		self.deltaUEigBuf)
	if self.fluxLimiter[0] > 0 then
		self.calcFluxKernel:setArg(6, self.rEigBuf)
	end

	-- TODO put this in solver/solver.lua ?
	if self.eqn.useSourceTerm then
		self.addSourceKernel = self.solverProgram:kernel'addSource'
		self.addSourceKernel:setArg(1, self.UBuf)
	end

	if self.checkFluxError or self.checkOrthoError then
		self.calcErrorsKernel = self.solverProgram:kernel(
			'calcErrors',
			self.errorBuf,
			self.waveBuf,
			self.eigenBuf)
	end
end

function Roe:addConvertToTexs()
	Roe.super.addConvertToTexs(self)

	for j,xj in ipairs(xNames) do
		self:addConvertToTex{
			name = 'wave',
			varCodePrefix = [[
	int indexInt = ]]..(j-1)..[[ + dim * index;
	const global real* wave = buf + indexInt * numWaves;
]],
			vars = range(0, self.eqn.numWaves-1):map(function(i)
				return {[xj..'_'..i] = '*value = wave['..i..'];'}
			end),
		}
	end

	local eigenDisplayVars = self.eqn:getEigenDisplayVars()
	if eigenDisplayVars and #eigenDisplayVars > 0 then
		for j,xj in ipairs(xNames) do
			self:addConvertToTex{
				name = 'eigen',
				type = self.eqn.eigen_t,
				varCodePrefix = [[
	int indexInt = ]]..(j-1)..[[ + dim * index;
	const global ]]..self.eqn.eigen_t..[[* eigen = buf + indexInt;
]],
				vars = table.map(eigenDisplayVars, function(kv)
					local k,v = next(kv)
					return {[xj..'_'..k] = v}
				end),
			}
		end
	end

	for j,xj in ipairs(xNames) do
		self:addConvertToTex{
			name = 'deltaUEig', 
			varCodePrefix = [[
	int indexInt = ]]..(j-1)..[[ + dim * index;
	const global real* deltaUEig = buf + indexInt * numWaves;
]],
			vars = range(0,self.eqn.numWaves-1):map(function(i)
				return {[xj..'_'..i] = '*value = deltaUEig['..i..'];'}
			end),
		}
	end

	if self.fluxLimiter[0] > 0 then
		for j,xj in ipairs(xNames) do
			self:addConvertToTex{
				name = 'rEig',
				varCodePrefix = [[
	int indexInt = ]]..(j-1)..[[ + dim * index;
	const global real* rEig = buf + indexInt * numWaves;
]],
				vars = range(0,self.eqn.numWaves-1):map(function(i)
					return {[xj..'_'..i] = '*value = rEig['..i..'];'}
				end),
			}
		end
	end
	
	-- TODO add kernels for each side
	if self.checkFluxError or self.checkOrthoError then	
		for j,xj in ipairs(xNames) do
			self:addConvertToTex{
				name = 'error', 
				varCodePrefix = [[
	int indexInt = ]]..(j-1)..[[ + dim * index;
]],
				useLog = true,
				type = 'error_t',
				vars = {
					{[xj..'_ortho'] = '*value = buf[indexInt].ortho;'},
					{[xj..'_flux'] = '*value = buf[indexInt].flux;'},
				},
			}
		end
	end
end

-- NOTICE this adds the contents of derivBuf and does not clear it
function Roe:calcDeriv(derivBuf, dt)
	self:boundary()
	
	if self.usePLM then
		self.calcLRKernel:setArg(2, ffi.new('real[1]', dt))
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcLRKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
	end

	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcEigenBasisKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}

	if self.checkFluxError or self.checkOrthoError then
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcErrorsKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
	end

	-- technically, if flux limiter isn't used, this can be merged into calcFlux (since no left/right reads need to be done)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDeltaUEigKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
	
	if self.fluxLimiter[0] > 0 then
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcREigKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
	end

	self.calcFluxKernel:setArg(5, ffi.new('real[1]', dt))
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcFluxKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}

	-- calcDerivFromFlux zeroes the derivative buffer
	self.calcDerivFromFluxKernel:setArg(0, derivBuf)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDerivFromFluxKernel, dim=self.dim, globalSize=self.globalSizeWithoutBorder:ptr(), localSize=self.localSize:ptr()}

	-- addSource adds to the derivative buffer
	if self.eqn.useSourceTerm then
		self.addSourceKernel:setArg(0, derivBuf)
		self.app.cmds:enqueueNDRangeKernel{kernel=self.addSourceKernel, dim=self.dim, globalSize=self.globalSizeWithoutBorder:ptr(), localSize=self.localSize:ptr()}
	end
end

return Roe

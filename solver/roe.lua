local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local math = require 'ext.math'
local file = require 'ext.file'
local template = require 'template'
local FiniteVolumeSolver = require 'solver.fvsolver'

local common = require 'common'()
local xNames = common.xNames
local symNames = common.symNames
local from3x3to6 = common.from3x3to6 
local from6to3x3 = common.from6to3x3 
local sym = common.sym


local Roe = class(FiniteVolumeSolver)
Roe.name = 'Roe'

function Roe:createBuffers()
	Roe.super.createBuffers(self)

	-- to get sizeof
	ffi.cdef(self.eqn:getEigenTypeCode())
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

	self.calcFluxKernelObj = self.solverProgramObj:kernel'calcFlux'
	self.calcFluxKernelObj.obj:setArg(1, self.fluxBuf)
end

local realptr = ffi.new'realparam[1]'
local function real(x)
	realptr[0] = x
	return realptr
end

-- NOTICE this adds the contents of derivBufObj and does not clear it
function Roe:calcDeriv(derivBufObj, dt)
if self.checkNaNs then assert(math.isfinite(dt)) end
	local dtArg = real(dt)

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
	
	if self.usePLM then
		self.calcLRKernelObj(self.solverBuf, self:getULRBuf(), self.UBuf, dtArg)
	end

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end

	self.calcFluxKernelObj.obj:setArg(0, self.solverBuf)
	self.calcFluxKernelObj.obj:setArg(2, self:getULRBuf())
	self.calcFluxKernelObj.obj:setArg(3, dtArg)
	self.calcFluxKernelObj()

if self.checkNaNs then assert(self:checkFinite(self.fluxBufObj)) end
if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end

-- [=[ this is from the 2017 Zingale book
	if self.useCTU then
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
		-- if we're using CTU then ...
		-- 1) calc fluxes based on a slope-limiter method (PLM, etc)
		-- 2) at each interface, integrate each dimension's LR states by all other dimensions' fluxes with a timestep of -dt/2
		--	( don't use the deriv buf because it already has the sum of all dimensions' flux differences)
		self.updateCTUKernelObj(self.solverBuf, self:getULRBuf(), self.fluxBuf, dtArg)
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end

		-- now we need to calcBounds on the ULR
		-- TODO this will break for mirror conditions
		-- because I haven't got the boundary code flexible enough to operate on specific fields within the L & R fields of the ULRBuf
		self:boundaryLR()

		-- 3) use the final LR states to calculate the flux ...

		-- the rest of this matches above
		-- maybe use 'repeat'?
		
--if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
--		self.calcEigenBasisKernelObj()
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
		self.calcFluxKernelObj()
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
if self.checkNaNs then assert(self:checkFinite(self.fluxBufObj)) end
	end
--]=]

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(self.fluxBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end
	
	self.calcDerivFromFluxKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivFromFluxKernelObj()

if self.checkNaNs then assert(self:checkFinite(self.UBufObj)) end
if self.checkNaNs then assert(self:checkFinite(derivBufObj)) end

end

return Roe

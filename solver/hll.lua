local ffi = require 'ffi'
local class = require 'ext.class'
local file = require 'ext.file'
local table = require 'ext.table'
local template = require 'template'
local FVSolver = require 'solver.fvsolver'

local HLL = class(FVSolver)
HLL.name = 'HLL'

HLL.solverCodeFile = 'solver/hll.cl'

--HLL.calcWaveMethod = 'Davis direct'
HLL.calcWaveMethod = 'Davis direct bounded'

function HLL:getSolverCode()
	return table{
		HLL.super.getSolverCode(self),
	
		-- before this went above solver/plm.cl, now it's going after it ...
		template(file[self.solverCodeFile], {
			solver = self,
			eqn = self.eqn,
			clnumber = require 'cl.obj.number',
		}),
	}:concat'\n'
end

function HLL:refreshSolverProgram()
	HLL.super.refreshSolverProgram(self)
	
	self.calcFluxKernelObj = self.solverProgramObj:kernel'calcFlux'
end

local realptr = ffi.new'realparam[1]'
local function real(x)
	realptr[0] = x
	return realptr
end
function HLL:calcDeriv(derivBufObj, dt)
	local dtArg = real(dt)
	
	if self.usePLM then
		self.calcLRKernelObj(self.solverBuf, self:getULRBuf(), self.UBuf, dtArg)
	end

	self.calcFluxKernelObj(self.solverBuf, self.fluxBuf, self:getULRBuf())

-- [=[ this is from the 2017 Zingale book
	if self.useCTU then	-- see solver/roe.lua for a description of why this is how this is
		self.updateCTUKernelObj(self.solverBuf, self:getULRBuf(), self.fluxBuf, dtArg)
	
		for _,obj in ipairs(self.lrBoundaryKernelObjs) do
			obj()
		end

		self.calcFluxKernelObj()
	end
--]=]
	
	self.calcDerivFromFluxKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcDerivFromFluxKernelObj()
end

return HLL

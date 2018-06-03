local ffi = require 'ffi'
local class = require 'ext.class'
local file = require 'ext.file'
local table = require 'ext.table'
local template = require 'template'
local FVSolver = require 'solver.fvsolver'

local HLL = class(FVSolver)
HLL.name = 'HLL'

function HLL:getSolverCode()
	return table{
		HLL.super.getSolverCode(self),
	
		-- before this went above solver/plm.cl, now it's going after it ...
		template(file['solver/hll.cl'], {
			solver = self,
			eqn = self.eqn,
			clnumber = require 'cl.obj.number',
		}),
	}:concat'\n'
end

function HLL:refreshSolverProgram()
	HLL.super.refreshSolverProgram(self)

	self.calcFluxKernelObj = self.solverProgramObj:kernel(
		'calcFlux',
		self.fluxBuf,
		self.getULRBuf)

	-- TODO put this in solver/gridsolver.lua ?
	if self.eqn.useSourceTerm then
		self.addSourceKernelObj = self.solverProgramObj:kernel{name='addSource', domain=self.domainWithoutBorder}
		self.addSourceKernelObj.obj:setArg(1, self.UBuf)
	end
end

local realptr = ffi.new'real[1]'
local function real(x)
	realptr[0] = x
	return realptr
end
function HLL:calcDeriv(derivBuf, dt)
	local dtArg = real(dt)
	
	self:boundary()
	
	if self.usePLM then
		self.calcLRKernelObj.obj:setArg(2, dtArg)
		self.calcLRKernelObj()
	end
	
	self.calcFluxKernelObj()

-- [=[ this is from the 2017 Zingale book
	if self.useCTU then	-- see solver/roe.lua for a description of why this is how this is
		self.updateCTUKernelObj.obj:setArg(2, dtArg)
		self.updateCTUKernelObj()
		
		self.lrBoundaryKernelObj()
		
		self.calcFluxKernelObj()
	end
--]=]
	
	self.calcDerivFromFluxKernelObj(derivBuf)
	
	-- addSource adds to the derivative buffer
	if self.eqn.useSourceTerm then
		-- can't use call because this uses the without-border size
		self.addSourceKernelObj(derivBuf)
	end
end

return HLL

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

	-- TODO put this in solver/solver.lua ?
	if self.eqn.useSourceTerm then
		self.addSourceKernelObj = self.solverProgramObj:kernel{name='addSource', domain=self.domainWithoutBorder}
		self.addSourceKernelObj.obj:setArg(1, self.UBuf)
	end
end

function HLL:calcDeriv(derivBuf, dt)
	self:boundary()
	
	if self.usePLM then
		self.calcLRKernelObj.obj:setArg(2, ffi.new('real[1]', dt))
		self.calcLRKernelObj()
	end
	
	self.calcFluxKernelObj()

	-- calcDerivFromFlux zeroes the derivative buffer
	self.calcDerivFromFluxKernelObj(derivBuf)

	-- addSource adds to the derivative buffer
	if self.eqn.useSourceTerm then
		-- can't use call because this uses the without-border size
		self.addSourceKernelObj(derivBuf)
	end
end

return HLL

local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local BSSNOKFiniteDifferenceEquation = require 'eqn.bssnok-fd'
local Solver = require 'solver.solver'

local BSSNOKFiniteDifferenceSolver = class(Solver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOKFiniteDifferenceSolver'

function BSSNOKFiniteDifferenceSolver:init(...)
	BSSNOKFiniteDifferenceSolver.super.init(self, ...)
	self.name = nil	-- don't append the eqn name to this
end

function BSSNOKFiniteDifferenceSolver:createEqn(eqn)
	self.eqn = BSSNOKFiniteDifferenceEquation(self)
end

function BSSNOKFiniteDifferenceSolver:refreshSolverProgram()
	BSSNOKFiniteDifferenceSolver.super.refreshSolverProgram(self)
	
	self.calcDerivKernel = self.solverProgram:kernel'calcDeriv'
	self.calcDerivKernel:setArg(1, self.UBuf)
end

function BSSNOKFiniteDifferenceSolver:refreshInitStateProgram()
	BSSNOKFiniteDifferenceSolver.super.refreshInitStateProgram(self)

	self.initConnUBarKernel = self.initStateProgram:kernel('init_connBarU', self.UBuf)
end

-- add an option for fixed Minkowsky boundary spacetime
function BSSNOKFiniteDifferenceSolver:createBoundaryOptions()
	BSSNOKFiniteDifferenceSolver.super.createBoundaryOptions(self)

	self.boundaryOptions:insert{
		fixed = function(args)
			local U = 'buf['..args.index'j'..']'
			return template([[
	<?=U?>.alpha = 1.;
	<?=U?>.beta_u = _real3(0,0,0);
	<?=U?>.gammaBar_ll = _sym3(1,0,0,1,0,1);
	<?=U?>.phi = 0;
	<?=U?>.K = 0;
	<?=U?>.ATilde_ll = _sym3(1,0,0,1,0,1);
	<?=U?>.connBar_u = _real3(0,0,0);
]], {U=U})
		end,
	}
end

function BSSNOKFiniteDifferenceSolver:resetState()
	BSSNOKFiniteDifferenceSolver.super.resetState(self)
	
	self.app.cmds:enqueueNDRangeKernel{kernel=self.initConnUBarKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
	self:boundary()
	self.app.cmds:finish()
end

function BSSNOKFiniteDifferenceSolver:getCalcDTCode() end
function BSSNOKFiniteDifferenceSolver:refreshCalcDTKernel() end
function BSSNOKFiniteDifferenceSolver:calcDT() return self.fixedDT end

function BSSNOKFiniteDifferenceSolver:calcDeriv(derivBuf, dt)
	self.calcDerivKernel:setArg(0, derivBuf)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDerivKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
end

return BSSNOKFiniteDifferenceSolver

local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local BSSNOKFiniteDifferenceEquation = require 'eqn.bssnok-fd'
local Solver = require 'solver.solver'

local BSSNOKFiniteDifferenceSolver = class(Solver)
BSSNOKFiniteDifferenceSolver.name = 'BSSNOKFiniteDifferenceSolver'

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

-- override and implement a constant boundary condition
-- TODO make the boundary conditions more flexible
-- right now there's no way to specify this in sell.app.boundaryMethods
function BSSNOKFiniteDifferenceSolver:refreshBoundaryProgram()
	self.boundaryProgram, self.boundaryKernel = 
		self:createBoundaryProgramAndKernel{
			type = self.eqn.cons_t,
			-- remap from enum/combobox int values to names
			methods = table.map(self.boundaryMethods, function(v,k)
				return function(U)
					return template([[
	<?=U?>.alpha = 1.;
	<?=U?>.beta_u = _real3(0,0,0);
	<?=U?>.gammaBar_ll = (sym3){.s={1,0,0,1,0,1}};
	<?=U?>.phi = 0;
	<?=U?>.K = 0;
	<?=U?>.ATilde_ll = (sym3){.s={1,0,0,1,0,1}};
	<?=U?>.connBar_u = _real3(0,0,0);
]],						{U = U})
				end, k
			end),
			mirrorVars = self.eqn.mirrorVars,
		}
	self.boundaryKernel:setArg(0, self.UBuf)
end


function BSSNOKFiniteDifferenceSolver:resetState()
	BSSNOKFiniteDifferenceSolver.super.resetState(self)
	
	self.app.cmds:enqueueNDRangeKernel{kernel=self.initConnUBarKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self:boundary()
	self.app.cmds:finish()
end

function BSSNOKFiniteDifferenceSolver:getCalcDTCode() end
function BSSNOKFiniteDifferenceSolver:refreshCalcDTKernel() end
function BSSNOKFiniteDifferenceSolver:calcDT()
	return self.fixedDT
end

function BSSNOKFiniteDifferenceSolver:calcDeriv(derivBuf, dt)
	self.calcDerivKernel:setArg(0, derivBuf)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDerivKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
end

return BSSNOKFiniteDifferenceSolver

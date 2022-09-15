local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local ig = require 'imgui'

-- TODO make this a ctor parameter 
local Poisson = require(
	cmdline.srhdSelfGravPoissonSolver 
	and 'hydro.op.poisson_'..cmdline.srhdSelfGravPoissonSolver
	or 'hydro.op.poisson_krylov'		-- Krylov
	--or 'hydro.op.poisson_jacobi'		-- Jacobi
)


local SRHDSelfGrav = class(Poisson)

SRHDSelfGrav.name = 'srhd_selfgrav'
SRHDSelfGrav.enableField = 'useGravity'

SRHDSelfGrav.guiVars = {
	{name='gravitationalConstant', value=2},
}

function SRHDSelfGrav:init(args)
	args.verbose = cmdline.selfGravVerbose
	args.linearSolver = cmdline.selfGravLinearSolver or 'conjres'	-- so far works best for selfgrav
	SRHDSelfGrav.super.init(self, args)
	self.solver[self.enableField] = not not self.solver[self.enableField]
end

function SRHDSelfGrav:getSymbolFields()
	return SRHDSelfGrav.super.getSymbolFields(self):append{
		'copyPotentialToReduce',	-- TODO use this?
		'offsetPotential',			-- TODO this too?
		'calcGravityAccel',
		'calcGravityDeriv',
	}
end

-- params for hydro/op/poisson.cl 
function SRHDSelfGrav:getPoissonDivCode()
	return self.solver.eqn:template([[
	source = 4. * M_PI * U->rho 
		* solver->gravitationalConstant / unit_m3_per_kg_s2;
]], {
		op = self,
	})
end

function SRHDSelfGrav:getPoissonCode()
	return file'hydro/op/srhd-selfgrav.cl':read()
end

function SRHDSelfGrav:initCodeModules()
	SRHDSelfGrav.super.initCodeModules(self)

	local solver = self.solver
	solver.solverModulesEnabled[self.symbols.calcGravityDeriv] = true
	--solver.solverModulesEnabled[self.symbols.copyPotentialToReduce] = true
	--solver.solverModulesEnabled[self.symbols.offsetPotential] = true
end

function SRHDSelfGrav:refreshSolverProgram()
	SRHDSelfGrav.super.refreshSolverProgram(self)
	
	local solver = self.solver
	self.calcGravityDerivKernelObj = solver.solverProgramObj:kernel(self.symbols.calcGravityDeriv)
	self.calcGravityDerivKernelObj.obj:setArg(0, solver.solverBuf)
	self.calcGravityDerivKernelObj.obj:setArg(2, solver.UBuf)
end

function SRHDSelfGrav:updateGUI()
	SRHDSelfGrav.super.updateGUI(self)
	ig.igPushID_Str'SRHDSelfGrav behavior'
	ig.luatableTooltipCheckbox('use gravity', self.solver, self.enableField)
	ig.igPopID()
end

function SRHDSelfGrav:addSource(derivBufObj)
	local solver = self.solver
	if not solver[self.enableField] then return end
		
	self:relax()
	self.calcGravityDerivKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcGravityDerivKernelObj()
end

return SRHDSelfGrav

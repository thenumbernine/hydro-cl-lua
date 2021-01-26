local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local table = require 'ext.table'
local class = require 'ext.class'
local tooltip = require 'hydro.tooltip'
local real = require 'hydro.real'	-- really 'realparam'

local half = require 'hydro.half'
local toreal, fromreal = half.toreal, half.fromreal


-- TODO make this a ctor parameter
local Poisson = require(
	cmdline.selfGravPoissonSolver 
	and 'hydro.op.poisson_'..cmdline.selfGravPoissonSolver
	-- Jacobi seems to be converging fastest atm. 
	-- however Jacobi seems to be most unstable for 3D.
	-- TODO multigrid or FFT?
	--or 'hydro.op.poisson_krylov'		-- Krylov. was working for a moment, but I broke it and now it is diverging.
	or 'hydro.op.poisson_jacobi'		-- Jacobi
)

-- TODO Schurr filter instead of 2n+1 point filter.

local SelfGrav = class(Poisson)

SelfGrav.name = 'selfgrav'
SelfGrav.enableField = 'useGravity'

-- source field we are using

function SelfGrav:init(args)
	-- TODO build super class based on what argument we chose?

	-- for krylov superclass:
	args.verbose = cmdline.selfGravVerbose
	args.linearSolver = cmdline.selfGravLinearSolver

	SelfGrav.super.init(self, args)
	
	self.solver[self.enableField] = not not self.solver[self.enableField]
end

function SelfGrav:getSymbolFields()
	return SelfGrav.super.getSymbolFields(self):append{
		'copyPotentialToReduce',
		'offsetPotential',
		'calcGravityAccel',
		'calcGravityDeriv',
	}
end

SelfGrav.guiVars = {
	{name='gravitationalConstant', value=1, units='m^3/(kg*s^2)'},
}

-- params for hydro/op/poisson.cl 
-- units of m^3/(kg*s^2) * kg/m^3 = 1/s^2
function SelfGrav:getPoissonDivCode()
	return self.solver.eqn:template([[
	source = 4. * M_PI * U->rho
		* solver->gravitationalConstant / unit_m3_per_kg_s2;	//'G'
]], {op=self})
end

function SelfGrav:getPoissonCode()
	return [[

//// MODULE_NAME: <?=calcGravityAccel?>

#define <?=calcGravityAccel?>(\
	/*real3 * const */accel_g,\
	/*constant <?=solver_t?> const * const */solver,\
	/*global <?=cons_t?> * const */U,\
	/*real3 const */pt\
) {\
	*(accel_g) = real3_zero;\
\
	<? for side=0,solver.dim-1 do ?>{\
		/* m/s^2 */\
		/* TODO grid coordinate influence? */\
		(accel_g)->s<?=side?> = (\
			U[solver->stepsize.s<?=side?>].<?=op.potentialField?> \
			- U[-solver->stepsize.s<?=side?>].<?=op.potentialField?>\
		) / (2. * solver->grid_dx.s<?=side?>);\
	}<? end ?>\
}

//// MODULE_NAME: <?=calcGravityDeriv?>
//// MODULE_DEPENDS: units realparam <?=calcGravityAccel?>

kernel void <?=calcGravityDeriv?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> * const derivBuffer,
	global <?=cons_t?> const * const UBuf,
	global <?=cell_t?> const * const cellBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	real3 const pt = cellBuf[index].pos;
	
	global <?=cons_t?> * const deriv = derivBuffer + index;
	global <?=cons_t?> const * const U = UBuf + index;

	real3 accel_g;
	<?=calcGravityAccel?>(&accel_g, solver, U, pt);

	// kg/(m^2 s) = kg/m^3 * m/s^2
	deriv->m = real3_sub(deriv->m, real3_real_mul(accel_g, U->rho));
	
	// kg/(m s^2) = (kg m^2 / s) * m/s^2
	deriv->ETotal -= real3_dot(U->m, accel_g);
}

//// MODULE_NAME: <?=copyPotentialToReduce?>

//TODO just use the display var kernels
kernel void <?=copyPotentialToReduce?>(
	constant <?=solver_t?> const * const solver,
	global real* reduceBuf,
	global const <?=cons_t?>* UBuf
) {
	<?=SETBOUNDS?>(0,0);
	reduceBuf[index] = UBuf[index].<?=op.potentialField?>;
}

//// MODULE_NAME: <?=offsetPotential?>

//keep potential energy negative
kernel void <?=offsetPotential?>(
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?>* UBuf,
	realparam ePotMax
) {
	<?=SETBOUNDS?>(0,0);
	global <?=cons_t?>* U = UBuf + index;
	U-><?=op.potentialField?> -= ePotMax;
}
]]
end

function SelfGrav:initCodeModules()
	SelfGrav.super.initCodeModules(self)
	
	local solver = self.solver
	solver.solverModulesEnabled[self.symbols.calcGravityDeriv] = true
	solver.solverModulesEnabled[self.symbols.copyPotentialToReduce] = true
	solver.solverModulesEnabled[self.symbols.offsetPotential] = true
end

function SelfGrav:refreshSolverProgram()
	SelfGrav.super.refreshSolverProgram(self)
	
	local solver = self.solver
	self.calcGravityDerivKernelObj = solver.solverProgramObj:kernel(self.symbols.calcGravityDeriv)
	self.calcGravityDerivKernelObj.obj:setArg(0, solver.solverBuf)
	self.calcGravityDerivKernelObj.obj:setArg(2, solver.UBuf)
	self.calcGravityDerivKernelObj.obj:setArg(3, solver.cellBuf)

	--TODO just use the display var kernels?
	self.copyPotentialToReduceKernelObj = solver.solverProgramObj:kernel(self.symbols.copyPotentialToReduce, solver.solverBuf, solver.reduceBuf, solver.UBuf)
	self.offsetPotentialKernelObj = solver.solverProgramObj:kernel(self.symbols.offsetPotential, solver.solverBuf, solver.UBuf)
end


-- TODO easier way to reduce min/max on display var ePot
--		but that means compiling the ePot dispaly var code even if the ePot display var is turned off ...
-- TODO stop calling them 'display var' since they're not just used for displaying
function SelfGrav:resetState()
	local solver = self.solver

	-- this does an initial relax()
	SelfGrav.super.resetState(self)

	if solver[self.enableField] then
		self:offsetPotential()
	end
end

function SelfGrav:updateGUI()
	SelfGrav.super.updateGUI(self)
	ig.igPushIDStr'SelfGrav behavior'
	tooltip.checkboxTable('use gravity', self.solver, self.enableField)
	ig.igPopID()
end

function SelfGrav:addSource(derivBufObj)
	local solver = self.solver
	if not solver[self.enableField] then return end

	self:relax()
	self.calcGravityDerivKernelObj.obj:setArg(1, derivBufObj.obj)
	self.calcGravityDerivKernelObj()

	-- while we are here, make sure the potential doesn't grow out of control
	self:offsetPotential()
end

function SelfGrav:offsetPotential()
	local solver = self.solver
	local ePotMin
	if self.verbose then
		self.copyPotentialToReduceKernelObj()
		ePotMin = fromreal(solver.reduceMin())
	end
	self.copyPotentialToReduceKernelObj()
	local ePotMax = fromreal(solver.reduceMax())
	
	self.offsetPotentialKernelObj.obj:setArg(2, real(ePotMax))
	self.offsetPotentialKernelObj()

	local new_ePotMin, new_ePotMax
	if self.verbose then
		self.copyPotentialToReduceKernelObj()
		new_ePotMin = fromreal(solver.reduceMin())
		self.copyPotentialToReduceKernelObj()
		new_ePotMax = fromreal(solver.reduceMax())
	
		print('offsetting potential energy from '..ePotMin..','..ePotMax..' to '..new_ePotMin..','..new_ePotMax)
	end
end

return SelfGrav

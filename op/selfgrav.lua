local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local tooltip = require 'tooltip'
local template = require 'template'

--local Poisson = require 'op.poisson'
local Poisson = require 'op.poisson_gmres'

local SelfGrav = class(Poisson)

SelfGrav.enableField = 'useGravity'

function SelfGrav:init(args)
	SelfGrav.super.init(self, args)
	self.solver[self.enableField] = not not self.solver[self.enableField]
end

SelfGrav.guiVars = {
	{name='gravitationalConstant', value=6.67384e-11},	-- m^3/(kg s^2)
}

-- params for op/poisson.cl 
-- units of m^3/(kg*s^2) * kg/m^3 = 1/s^2
function SelfGrav:getPoissonDivCode()
	return template([[
	source = solver->gravitationalConstant * unit_m3_per_kg_s2 * U->rho;
]], {
		op = self,
	})
end

function SelfGrav:getPoissonCode()
	return template([[
kernel void calcGravityDeriv(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* derivBuffer,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	
	global <?=eqn.cons_t?>* deriv = derivBuffer + index;
	const global <?=eqn.cons_t?>* U = UBuf + index;
		
	//for (int side = 0; side < dim; ++side) {
	<? for side=0,solver.dim-1 do ?>{
		const int side = <?=side?>;
		int indexL = index - solver->stepsize.s<?=side?>;
		int indexR = index + solver->stepsize.s<?=side?>;

		// m/s^2
		real accel_g = (
			UBuf[indexR].<?=op.potentialField?> 
			- UBuf[indexL].<?=op.potentialField?>
		) / (2. * cell_dx<?=side?>(x));

		// kg/(m^2 s) = kg/m^3 * m/s^2
		deriv->m.s[side] -= U->rho * accel_g;
	
		// kg/(m s^2) = (kg m^2 / s) * m/s^2
		deriv->ETotal -= U->m.s[side] * accel_g;
	}<? end ?>
}

//TODO just use the display var kernels
kernel void copyPotentialToReduce(
	constant <?=solver.solver_t?>* solver,
	global real* reduceBuf,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	reduceBuf[index] = UBuf[index].<?=op.potentialField?>;
}

//hmm, if ePot is negative then we get the cool turbulence effect
//but if it is positive, esp >1, then we get no extra behavior ... a static sphere
//and if it is too big then it explodes
kernel void offsetPotentialAndAddToTotal(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	realparam ePotMin
) {
	const real basePotential = 0.;
	//const real basePotential = 1.;
	
	SETBOUNDS(0,0);
	global <?=eqn.cons_t?>* U = UBuf + index;
	U-><?=op.potentialField?> += basePotential - ePotMin;
	real source = 0.;
<?=op:getPoissonDivCode()?>
	U->ETotal += source * U-><?=op.potentialField?>;
}
]], {
		op = self,
		solver = self.solver,
		eqn = self.solver.eqn,
	})
end

function SelfGrav:refreshSolverProgram()
	SelfGrav.super.refreshSolverProgram(self)
	
	local solver = self.solver
	self.calcGravityDerivKernelObj = solver.solverProgramObj:kernel'calcGravityDeriv'
	self.calcGravityDerivKernelObj.obj:setArg(0, solver.solverBuf)
	self.calcGravityDerivKernelObj.obj:setArg(2, solver.UBuf)

	--TODO just use the display var kernels?
	self.copyPotentialToReduceKernelObj = solver.solverProgramObj:kernel('copyPotentialToReduce', solver.solverBuf, solver.reduceBuf, solver.UBuf)
	self.offsetPotentialAndAddToTotalKernelObj = solver.solverProgramObj:kernel('offsetPotentialAndAddToTotal', solver.solverBuf, solver.UBuf)
end

local realptr = ffi.new'realparam[1]'
local function real(x)
	realptr[0] = x
	return realptr
end


-- TODO easier way to reduce min/max on display var ePot
--		but that means compiling the ePot dispaly var code even if the ePot display var is turned off ...
-- TODO stop calling them 'display var' since they're not just used for displaying
function SelfGrav:resetState()
	local solver = self.solver

	-- this does an initial relax()
	SelfGrav.super.resetState(self)

	self.copyPotentialToReduceKernelObj()

	local ePotMin = solver.reduceMin()
	self.copyPotentialToReduceKernelObj()
	local ePotMax = solver.reduceMax()

	self.offsetPotentialAndAddToTotalKernelObj.obj:setArg(2, real(ePotMin))
	self.offsetPotentialAndAddToTotalKernelObj()
	
	self.copyPotentialToReduceKernelObj()
	local new_ePotMin = solver.reduceMin()
	self.copyPotentialToReduceKernelObj()
	local new_ePotMax = solver.reduceMax()

	print('offsetting potential energy from '..ePotMin..','..ePotMax..' to '..new_ePotMin..','..new_ePotMax)
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
end

return SelfGrav

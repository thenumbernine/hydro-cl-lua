local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local tooltip = require 'tooltip'
local template = require 'template'

local Poisson = require 'op.poisson'
--local Poisson = require 'op.poisson_gmres'

local SelfGrav = class(Poisson)

-- potential field is inhertied from Poisson: ePot
SelfGrav.densityField = 'rho'

SelfGrav.enableField = 'useGravity'

SelfGrav.gravitationConstant = 1	---- 6.67384e-11 m^3 / (kg s^2)

function SelfGrav:init(args)
	SelfGrav.super.init(self, args)
	self.densityField = args.densityField	
	self.solver[self.enableField] = not not self.solver[self.enableField]
end

-- params for op/poisson.cl 
function SelfGrav:getCalcRhoCode()
	return template([[
	//maybe a 4pi?  or is that only in the continuous case?
	rho = <?=clnumber(op.gravitationConstant)?>
		* U-><?=op.densityField?>;
]], {
		op = self,
		solver = self.solver,
		eqn = self.solver.eqn,
		clnumber = require 'cl.obj.number',
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
		int indexL = index - stepsize.s<?=side?>;
		int indexR = index + stepsize.s<?=side?>;
	
		real gravity = (UBuf[indexR].<?=self.potentialField?> - UBuf[indexL].<?=self.potentialField?>) / (2. * dx<?=side?>_at(i));

		deriv->m.s[side] -= U-><?=self.densityField?> * gravity;
		deriv->ETotal -= U-><?=self.densityField?> * gravity * U->m.s[side];
	}<? end ?>
}

//TODO just use the display var kernels
kernel void reduce_ePot(
	global real* reduceBuf,
	global const <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	reduceBuf[index] = UBuf[index].<?=self.potentialField?>;
}

//hmm, if ePot is negative then we get the cool turbulence effect
//but if it is positive, esp >1, then we get no extra behavior ... a static sphere
//and if it is too big then it explodes
kernel void offsetPotentialAndAddToTotal(
	global <?=eqn.cons_t?>* UBuf,
	realparam ePotMin
) {
	const real basePotential = 0.;
	//const real basePotential = 1.;
	
	SETBOUNDS(0,0);
	global <?=eqn.cons_t?>* U = UBuf + index;
	U-><?=self.potentialField?> += basePotential - ePotMin;
	U->ETotal += U-><?=self.densityField?> * U-><?=self.potentialField?>;
}
]], {
		self = self,
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

	--TODO just use the display var kernels
	self.reduce_ePotKernelObj = solver.solverProgramObj:kernel('reduce_ePot', solver.reduceBuf, solver.UBuf)
	
	self.offsetPotentialAndAddToTotalKernelObj = solver.solverProgramObj:kernel('offsetPotentialAndAddToTotal', solver.UBuf)
end

local realptr = ffi.new'realparam[1]'
local function real(x)
	realptr[0] = x
	return realptr
end

function SelfGrav:resetState()
	SelfGrav.super.resetState(self)

	local solver = self.solver
	
	-- TODO easier way to reduce min/max on display var ePot
	--		but that means compiling the ePot dispaly var code even if the ePot display var is turned off ...
	-- TODO stop calling them 'display var' since they're not just used for displaying
	self.reduce_ePotKernelObj()
	local ePotMin = solver.reduceMin()
	self.reduce_ePotKernelObj()
	local ePotMax = solver.reduceMax()

	self.offsetPotentialAndAddToTotalKernelObj.obj:setArg(1, real(ePotMin))
	self.offsetPotentialAndAddToTotalKernelObj()
	
	self.reduce_ePotKernelObj()
	local new_ePotMin = solver.reduceMin()
	self.reduce_ePotKernelObj()
	local new_ePotMax = solver.reduceMax()

print('offsetting potential energy from '..ePotMin..','..ePotMax..' to '..new_ePotMin..','..new_ePotMax)
end

function SelfGrav:updateGUI()
	SelfGrav.super.updateGUI(self)
	ig.igPushIDStr'SelfGrav behavior'
	tooltip.checkboxTable('use gravity', self.solver, self.enableField)
	ig.igPopID()
end

function SelfGrav:step(dt)
	local solver = self.solver
	if not solver[self.enableField] then return end
	solver.integrator:integrate(dt, function(derivBuf)
		self:relax()
		self.calcGravityDerivKernelObj.obj:setArg(1, derivBuf)
		self.calcGravityDerivKernelObj()
	end)
end

return SelfGrav

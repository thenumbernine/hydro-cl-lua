local class = require 'ext.class'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'

local NavierStokesDivFree = class(Equation)
NavierStokesDivFree.name = 'Navier-Stokes div-free' 
NavierStokesDivFree.hasCalcDTCode = true
NavierStokesDivFree.numStates = 4

NavierStokesDivFree.reflectVars = {
	mirror = {
		{'v.x'},
		{'v.y'},
		{'v.z'},
	},
}

NavierStokesDivFree.initStates = require 'init.euler'

function NavierStokesDivFree:getTypeCode()
	return template([[
typedef union {
	real ptr[4];
	struct {
		real rho;
		real3 v;
	};
} <?=eqn.prim_t?>;

typedef <?=eqn.prim_t?> <?=eqn.cons_t?>; 
]], {
		eqn = self,
	})
end

NavierStokesDivFree.initStateCode = [[
kernel void initState(
	global <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_real_mul(real3_add(mins, maxs), .5);
	bool lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
	real rho = 0;
	real3 v = real3_zero;
	
	// throw-away:
	real P = 0;
	real3 B = real3_zero;
	real ePot = 0;

	<?=code?>

	<?=eqn.cons_t?> U = {
		.rho = rho,
		.v = v,
	};
	UBuf[index] = U;
}
]]

NavierStokesDivFree.solverCodeFile = 'eqn/euler.cl'

function NavierStokesDivFree:getDisplayVars()
	local vars = table{
		{rho = 'value.vreal = U->rho;'},
		{v = 'value.vreal3 = U->v;', type='real3'},
		{m = 'value = real3_real_mul(U->v, U->rho);', type='real3'},
		--{P = 'value.vreal = W.P;'},
		--{eInt = 'value.vreal = calc_eInt(W);'},
		--{eKin = 'value.vreal = calc_eKin(W);'},
		--{ePot = 'value.vreal = U->ePot;'},
		--{eTotal = 'value.vreal = U->ETotal / W.rho;'},
		--{EInt = 'value.vreal = calc_EInt(W);'},
		{EKin = 'value.vreal = .5 * U->rho * coordLen(U->v, x);'},
		--{EPot = 'value.vreal = U->rho * U->ePot;'},
		--{ETotal = 'value.vreal = U->ETotal;'},
		--{S = 'value.vreal = W.P / pow(W.rho, (real)heatCapacityRatio);'},
		--{H = 'value.vreal = calc_H(W.P);'},
		--{h = 'value.vreal = calc_h(W.rho, W.P);'},
		--{HTotal = 'value.vreal = calc_HTotal(W.P, U->ETotal);'},
		--{hTotal = 'value.vreal = calc_hTotal(W.rho, W.P, U->ETotal);'},
		--{['Speed of Sound'] = 'value.vreal = calc_Cs(&W);'},
		--{['Mach number'] = 'value.vreal = coordLen(W.v, x) / calc_Cs(&W);'},
	}
	
	vars:insert(self:createDivDisplayVar{
		field = 'v', 
		units = 'kg/(m^3*s)',
	})

	vars:insert(self:createCurlDisplayVar{
		field = 'v',
		units = 'm/s^2',
	})

	return vars
end

return NavierStokesDivFree

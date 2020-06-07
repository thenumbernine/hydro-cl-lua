local class = require 'ext.class'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'hydro.eqn.eqn'

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

NavierStokesDivFree.initStates = require 'hydro.init.euler'

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

NavierStokesDivFree.solverCodeFile = 'hydro/eqn/euler.cl'

function NavierStokesDivFree:getDisplayVars()
	local vars = table{
		{name='rho', code='value.vreal = U->rho;'},
		{name='v', code='value.vreal3 = U->v;', type='real3'},
		{name='m', code='value = real3_real_mul(U->v, U->rho);', type='real3'},
		--{name='P', code='value.vreal = W.P;'},
		--{name='eInt', code='value.vreal = calc_eInt(W);'},
		--{name='eKin', code='value.vreal = calc_eKin(W);'},
		--{name='ePot', code='value.vreal = U->ePot;'},
		--{name='eTotal', code='value.vreal = U->ETotal / W.rho;'},
		--{name='EInt', code='value.vreal = calc_EInt(W);'},
		{name='EKin', code='value.vreal = .5 * U->rho * coordLen(U->v, x);'},
		--{name='EPot', code='value.vreal = U->rho * U->ePot;'},
		--{name='ETotal', code='value.vreal = U->ETotal;'},
		--{name='S', code='value.vreal = W.P / pow(W.rho, (real)heatCapacityRatio);'},
		--{name='H', code='value.vreal = calc_H(W.P);'},
		--{name='h', code='value.vreal = calc_h(W.rho, W.P);'},
		--{name='HTotal', code='value.vreal = calc_HTotal(W.P, U->ETotal);'},
		--{name='hTotal', code='value.vreal = calc_hTotal(W.rho, W.P, U->ETotal);'},
		--{name='Speed of Sound' code='value.vreal = calc_Cs(&W);'},
		--{name='Mach number' code='value.vreal = coordLen(W.v, x) / calc_Cs(&W);'},
	}
	
	vars:insert(self:createDivDisplayVar{field='v', units='1/s'} or nil)
	vars:insert(self:createCurlDisplayVar{field='v', units='1/s'} or nil)

	return vars
end

return NavierStokesDivFree
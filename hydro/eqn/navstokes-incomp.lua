local class = require 'ext.class'
local Equation = require 'hydro.eqn.eqn'

local NavierStokesDivFree = class(Equation)
NavierStokesDivFree.name = 'Navier-Stokes div-free' 

NavierStokesDivFree.numStates = 4

NavierStokesDivFree.initConds = require 'hydro.init.euler':getList()

NavierStokesDivFree.consVars = table{
	{name='rho', type='real', units='kg/m^3'},
	{name='v', type='real3', units='m/s', variance='u'},
}

-- hmm, if this is using euler.cl then how can we extend this for applyInitCondCell?
-- I might have to bring the getModuleCode_applyInitCondCell function back 
NavierStokesDivFree.initCondCode = [[
void applyInitCondCell(
	constant <?=solver_t?> const * const solver,
	constant <?=initCond_t?> const * const initCond,
	global <?=cons_t?> * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;
	real3 const mids = real3_real_mul(real3_add(mins, maxs), .5);
	bool const lhs = x.x < mids.x
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

	*U = (<?=cons_t?>){
		.rho = rho,
		.v = v,
	};
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

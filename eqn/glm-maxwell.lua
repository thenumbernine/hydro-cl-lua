--[[ 2000 Munz
But then I changed the eigendecomposition to match my Trangenstein one with epsilons and mus
And then I changed epsilon E -> D, so it's based on B and D 
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local Equation = require 'eqn.eqn'
local clnumber = require 'cl.obj.number'
local template = require 'template'

local GLM_Maxwell = class(Equation)
GLM_Maxwell.name = 'GLM_Maxwell'
GLM_Maxwell.numWaves = 8
GLM_Maxwell.numIntStates = 8

GLM_Maxwell.consVars = {
	{D = 'real3'},
	{B = 'real3'},
	{phi = 'real'},
	{psi = 'real'},
	{rhoCharge = 'real'},
	{sigma = 'real'},
	{sqrt_1_eps = 'real'},
	{sqrt_1_mu = 'real'},
}

GLM_Maxwell.mirrorVars = {{'D.x', 'B.x'}, {'D.y', 'B.y'}, {'D.z', 'B.z'}}

GLM_Maxwell.hasEigenCode = true
GLM_Maxwell.hasFluxFromConsCode = true
GLM_Maxwell.useSourceTerm = true
GLM_Maxwell.roeUseFluxFromCons = true

GLM_Maxwell.initStates = require 'init.euler'

GLM_Maxwell.postComputeFluxCode = template([[
		//TODO shouldn't I be transforming both the left and right fluxes by the metrics at their respective coordinates?
		//flux is computed raised via Levi-Civita upper
		//so here we lower it
		real _1_sqrt_det_g = 1. / sqrt_det_g_grid(x);
		flux.D = real3_real_mul(coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = real3_real_mul(coord_lower(flux.B, x), _1_sqrt_det_g);
]], {eqn=Maxwell})


function GLM_Maxwell:getCommonFuncCode()
	return template([[
real3 calc_E(<?=eqn.cons_t?> U) { 
	return real3_real_mul(U.D, U.sqrt_1_eps * U.sqrt_1_eps);
}
real3 calc_H(<?=eqn.cons_t?> U) { 
	return real3_real_mul(U.B, U.sqrt_1_mu * U.sqrt_1_mu);
}

real ESq(<?=eqn.cons_t?> U, real3 x) { return coordLenSq(calc_E(U), x); }
real BSq(<?=eqn.cons_t?> U, real3 x) { return coordLenSq(U.B, x); }
]], {
	eqn = self,
})
end

GLM_Maxwell.initStateCode = [[
kernel void initState(
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
	global <?=eqn.cons_t?>* U = UBuf + index;

	//used
	real3 E = real3_zero;
	real3 B = real3_zero;
	real conductivity = 1.;
	
	real permittivity = 1.;
	real permeability = 1.;
	
	//throw-away
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real ePot = 0;
	
	<?=code?>
	
	U->D = cartesianToCoord(real3_real_mul(E, permittivity), x);
	U->B = cartesianToCoord(B, x);
	U->phi = 0;
	U->psi = 0;
	U->sigma = conductivity;
	U->rhoCharge = 0;
	U->sqrt_1_eps = sqrt(1. / permittivity);
	U->sqrt_1_mu = sqrt(1. / permeability);
}
]]

GLM_Maxwell.solverCodeFile = 'eqn/glm-maxwell.cl'

-- k is 0,1,2
local function curl(eqn,k,result,field)
	local xs = {'x','y','z'}
	local i = (k+1)%3
	local j = (i+1)%3
	return {['curl '..field..' '..xs[k+1]] = template([[
	if (OOB(1,1)) {
		<?=result?> = 0.;
	} else {

<? if i+1 <= solver.dim then ?>
		global const <?=eqn.cons_t?>* Uim = U - stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Uip = U + stepsize.s<?=i?>;
		real vim_j = Uim-><?=field?>.s<?=j?>;
		real vip_j = Uip-><?=field?>.s<?=j?>;
<? else ?>
		real vim_j = 0.;
		real vip_j = 0.;
<? end?>

<? if j+1 <= solver.dim then ?>
		global const <?=eqn.cons_t?>* Ujm = U - stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + stepsize.s<?=j?>;
		real vjm_i = Ujm-><?=field?>.s<?=i?>;
		real vjp_i = Ujp-><?=field?>.s<?=i?>;
<? else ?>
		real vjm_i = 0.;
		real vjp_i = 0.;
<? end ?>

		<?=result?> = (vjp_i - vjm_i) / (2. * grid_dx<?=i?>)
				- (vip_j - vim_j) / (2. * grid_dx<?=j?>);
	}
]], {
		i = i,
		j = j,
		eqn = eqn,
		solver = eqn.solver,
		result = result,
		field = field,
	})}
end

function GLM_Maxwell:getDisplayVars()
	local vars = GLM_Maxwell.super.getDisplayVars(self):append{ 
		{S = [[
	*value_real3 = real3_cross(calc_E(*U), calc_H(*U));
]], type='real3'},
		{energy = [[
	*value = .5 * (coordLenSq(U->D, x) + coordLenSq(calc_H(*U), x));
]]},
	}:append(table{'D','B'}:map(function(var,i)
		local field = assert( ({D='D', B='B'})[var] )
		return {['div '..var] = template([[
	*value = .5 * (0.
<?
for j=0,solver.dim-1 do
?>		+ (U[stepsize.s<?=j?>].<?=field?>.s<?=j?> 
			- U[-stepsize.s<?=j?>].<?=field?>.s<?=j?>
		) / grid_dx<?=j?>
<?
end 
?>	);
]], {solver=self.solver, field=field})}
	end))

	for _,field in ipairs{'D', 'B'} do
		local v = range(0,2):map(function(i) 
			return curl(self,i,'value_real3->s'..i,field) 
		end)
		vars:insert{['curl '..field]= template([[
	<? for i=0,2 do ?>{
		<?=select(2,next(v[i+1]))?>
	}<? end ?>
]], {v=v}), type='real3'}
	end

	return vars
end

GLM_Maxwell.eigenVars = table{
	{sqrt_1_eps = 'real'},
	{sqrt_1_mu = 'real'},
}

function GLM_Maxwell:eigenWaveCodePrefix(side, eig, x, waveIndex)
	return template([[
	real v_p = <?=eig?>.sqrt_1_eps * <?=eig?>.sqrt_1_mu;
]], {
		eqn = self,
		eig = '('..eig..')',
	})
end


function GLM_Maxwell:eigenWaveCode(side, eig, x, waveIndex)
	return ({
		'-v_p * divPhiWavespeed',
		'-v_p * divPsiWavespeed',
		'-v_p',
		'-v_p',
		'v_p',
		'v_p',
		'v_p * divPhiWavespeed',
		'v_p * divPsiWavespeed',
	})[waveIndex+1] or error('got a bad waveIndex: '..waveIndex)
end

function GLM_Maxwell:consWaveCodePrefix(side, U, x, waveIndex) 
	return template([[
	real v_p = <?=U?>.sqrt_1_eps * <?=U?>.sqrt_1_mu;
]], {
		eqn = self,
		U = '('..U..')',
	})
end
GLM_Maxwell.consWaveCode = GLM_Maxwell.eigenWaveCode

return GLM_Maxwell

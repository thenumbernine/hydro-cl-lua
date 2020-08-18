--[[
based on Trangenstein

for curved space, I'll keep my vectors in covariant form
this way the Levi-Civita tensor in teh flux, multiplied by the connection coefficents, when used in a holonomic basis, makes them cancel
and when used in a static grid, the g^ij_,t terms vanish,
and you get

eps0 E_i,t - 1/sqrt(g) g_il 1/mu0  epsBar^ljk      B_k,j = -j_i
     B_i,t + 1/sqrt(g) g_il 1/eps0 epsBar^ljk eps0 E_k,j =   0

so during flux computation, I need to not apply the sqrt det g
and after computing the flux vector I need to apply the metric

How to add in the D and H fields...

I am really tempted to change eps0 E -> E
so that I can more freely mess with the aux fields:

in fact, in materials, D is a better candidate anyways, since formula are in D,t and B,t, and D = epsilon E, so using eps0 E is a good start

D_i,t - 1/sqrt(g) g_il epsBar^ljk  1/mu B_k,j = 1/sqrt(g) g_il epsBar^ljk  (1/mu)_,j B_k - J_i
B_i,t + 1/sqrt(g) g_il epsBar^ljk 1/eps D_k,j = 1/sqrt(g) g_il epsBar^ljk (1/eps)_,j B_k

D_i,t - 1/sqrt(g) g_il epsBar^ljk  1/mu B_k,j = -1/sqrt(g) g_il epsBar^ljk B_j  (1/mu)_,k - J_i
B_i,t + 1/sqrt(g) g_il epsBar^ljk 1/eps D_k,j = -1/sqrt(g) g_il epsBar^ljk D_j (1/eps)_,k

TODO now I need to add source terms of the permittivity and permeability gradients ...
that will look like ...

D_i,t - 1/sqrt(g) g_il epsBar^ljk  (1/mu)_k^l B_l,j = 1/sqrt(g) g_il epsBar^ljk  (1/mu)_k^l_,j B_l - J_i
B_i,t + 1/sqrt(g) g_il epsBar^ljk (1/eps)_k^l D_l,j = 1/sqrt(g) g_il epsBar^ljk (1/eps)_k^l_,j B_l

--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'

local Maxwell = class(Equation)
Maxwell.name = 'Maxwell'

Maxwell.roeUseFluxFromCons = true
Maxwell.useSourceTerm = true

-- don't incorporate the Conn^k_ij E_k terms into the flux
Maxwell.weightFluxByGridVolume = false

Maxwell.initConds = require 'hydro.init.euler'

--[[
args:
	scalar = 'real' or 'cplx'
--]]
function Maxwell:init(args)
	self.scalar = (args and args.scalar) or 'real'
	
	self.vec3 = self.scalar..'3'
	self.mat3x3 = self.scalar..'3x3'

	-- TODO tensor susceptibilty support ... but that affects the eigendecomposition ...
	self.susc_t = self.scalar
	--self.susc_t = self.mat3x3

	self.numRealsInScalar = ffi.sizeof(self.scalar) / ffi.sizeof'real'

	self.numIntStates = 6 * self.numRealsInScalar
	self.numWaves = 6 * self.numRealsInScalar

	self.consVars = {
		{name='D', type=self.vec3, units='C/m^2', variance='l'},		-- D_i
		{name='B', type=self.vec3, units='kg/(C*s)', variance='l'},		-- B_i
		{name='phi', type=self.scalar, units='C/m^2'},
		{name='psi', type=self.scalar, units='kg/(C*s)'},
		{name='rhoCharge', type=self.scalar, units='C/m^3'},
		{name='J', type=self.vec3, units='C/(m^2*s)'},
		{name='sqrt_1_eps', type=self.susc_t, units='(kg*m^3)^.5/(C*s)'},
		{name='sqrt_1_mu', type=self.susc_t, units='C/(kg*m)^.5'},
	}
	
	self.eigenVars = table{
		{name='sqrt_1_eps', type=self.susc_t, units='(kg*m^3)^.5/(C*s)'},
		{name='sqrt_1_mu', type=self.susc_t, units='C/(kg*m)^.5'},
	}

	Maxwell.super.init(self, args)

-- [=[ TODO combine this into the flux and remove this variable from calcDerivFV
	self.postComputeFluxCode = self:template[[
<? local vec3 = eqn.vec3 ?>
		//TODO shouldn't I be transforming both the left and right fluxes by the metrics at their respective coordinates?
		//flux is computed raised via Levi-Civita upper
		//so here we lower it
		real _1_sqrt_det_g = 1. / coord_sqrt_det_g(x);
		flux.D = <?=vec3?>_real_mul(eqn_coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = <?=vec3?>_real_mul(eqn_coord_lower(flux.B, x), _1_sqrt_det_g);
]]
--]=]

	local NoDiv = require 'hydro.op.nodiv'{
		poissonSolver = require 'hydro.op.poisson_jacobi',
	}

	self.solver.ops:insert(NoDiv{
		solver = self.solver,
		scalar = self.scalar,
		vectorField = 'B',
		potentialField = 'psi',
	})

	self.solver.ops:insert(NoDiv{
		solver = self.solver,
		scalar = self.scalar,
		vectorField = 'D',
		potentialField = 'phi',
		chargeField = 'rhoCharge',
	})
end

function Maxwell:initCodeModule_fluxFromCons()
	self.solver.modules:add{
		name = 'fluxFromCons',
		depends = {
			'solver.solver_t',
			'coord.normal',
			'eqn.cons_t',
			'eqn.prim_t',
			'eqn.common',	-- calc_E, calc_H
		},
		code = self:template[[
<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normalInfo_t n
) {
	<?=vec3?> E = calc_E(U);
	<?=vec3?> H = calc_H(U);
	<?=eqn.cons_t?> F;
	if (n.side == 0) {
		F.D = _<?=vec3?>(<?=zero?>, H.z, <?=neg?>(H.y));
		F.B = _<?=vec3?>(<?=zero?>, <?=neg?>(E.z), E.y);
	} else if (n.side == 1) {
		F.D = _<?=vec3?>(<?=neg?>(H.z), <?=zero?>, H.x);
		F.B = _<?=vec3?>(E.z, <?=zero?>, <?=neg?>(E.x));
	} else if (n.side == 2) {
		F.D = _<?=vec3?>(H.y, <?=neg?>(H.x), <?=zero?>);
		F.B = _<?=vec3?>(<?=neg?>(E.y), E.x, <?=zero?>);
	}
	F.phi = <?=zero?>;
	F.psi = <?=zero?>;
	F.D = <?=vec3?>_zero;
	F.rhoCharge = <?=zero?>;
	F.sqrt_1_eps = <?=susc_t?>_zero;
	F.sqrt_1_mu = <?=susc_t?>_zero;
	return F;
}
]],
	}
end

function Maxwell:initCodeModuleCommon()
	self.solver.modules:add{
		name = 'eqn.common',
		code = self:template[[
<? if scalar == 'real' then ?>

#define eqn_coordLenSq coordLenSq
#define eqn_cartesianToCoord cartesianToCoord
#define eqn_coord_lower coord_lower

<? elseif scalar == 'cplx' then ?>

real eqn_coordLenSq(cplx3 v, real3 x) {
	return coordLenSq(cplx3_re(v), x)
		+ coordLenSq(cplx3_im(v), x);
}

cplx3 eqn_cartesianToCoord(cplx3 v, real3 x) {
	return cplx3_from_real3_real3(
		cartesianToCoord(cplx3_re(v), x),
		cartesianToCoord(cplx3_im(v), x));
}

cplx3 eqn_coord_lower(cplx3 v, real3 x) {
	return cplx3_from_real3_real3(
		coord_lower(cplx3_re(v), x),
		coord_lower(cplx3_im(v), x));
}

<? end -- scalar ?>

<?=vec3?> calc_E(<?=eqn.cons_t?> U) { 
	return <?=vec3?>_<?=susc_t?>_mul(U.D, <?=susc_t?>_mul(U.sqrt_1_eps, U.sqrt_1_eps));
}
<?=vec3?> calc_H(<?=eqn.cons_t?> U) { 
	return <?=vec3?>_<?=susc_t?>_mul(U.B, <?=susc_t?>_mul(U.sqrt_1_mu, U.sqrt_1_mu));
}
]],
	}
end

Maxwell.initCondCode = [[
<? 
local cons_t = eqn.cons_t
local susc_t = eqn.susc_t
local scalar = eqn.scalar
local vec3 = eqn.vec3
local zero = scalar..'_zero'
?>

kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
	global <?=cons_t?>* U = UBuf + index;

	//used
	<?=vec3?> D = <?=vec3?>_zero;
	<?=vec3?> B = <?=vec3?>_zero;
	<?=scalar?> conductivity = <?=scalar?>_from_real(1.);
	<?=susc_t?> permittivity = <?=susc_t?>_from_real(1.);
	<?=susc_t?> permeability = <?=susc_t?>_from_real(1.);
	<?=scalar?> rhoCharge = <?=zero?>;

	//throw-away
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real ePot = 0;
	
<?=code?>
	
	U->D = eqn_cartesianToCoord(D, x);
	U->B = eqn_cartesianToCoord(B, x);
	U->phi = <?=zero?>;
	U->psi = <?=zero?>;
	U->J = <?=vec3?>_zero;
	U->rhoCharge = rhoCharge;
	U->sqrt_1_eps = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permittivity));
	U->sqrt_1_mu = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permeability));
}
]]

function Maxwell:getModuleDependsApplyInitCond()
	return {'eqn.common'}
end

function Maxwell:getModuleDependsSolver()
	return {
		'eqn.common',
		'coord_lower',
		'fluxFromCons',
	}
end

Maxwell.solverCodeFile = 'hydro/eqn/maxwell.cl'

function Maxwell:getEnv()
	local scalar = self.scalar
	local env = Maxwell.super.getEnv(self)
	env.vec3 = self.vec3
	env.susc_t = self.susc_t
	env.scalar = scalar
	env.zero = scalar..'_zero'
	env.inv = scalar..'_inv'
	env.neg = scalar..'_neg'
	env.fromreal = scalar..'_from_real'
	env.add = scalar..'_add'
	env.sub = scalar..'_sub'
	env.mul = scalar..'_mul'
	env.mul3 = scalar..'_mul3'
	env.real_mul = scalar..'_real_mul'
	env.sqrt = scalar..'_sqrt'
	env.abs = scalar..'_abs'
	return env
end

Maxwell.predefinedDisplayVars = {
	'U D',
	'U div D',
	'U phi',
	'U B',
	'U div B',
	'U psi',
}

function Maxwell:getDisplayVars()
	local env = self:getEnv()
	local vars = Maxwell.super.getDisplayVars(self)
	vars:append{ 
		{name='E', code=self:template[[	value.v<?=vec3?> = calc_E(*U);]], type=env.vec3, units='(kg*m)/(C*s)'},
		{name='H', code=self:template[[	value.v<?=vec3?> = calc_H(*U);]], type=env.vec3, units='C/(m*s)'},
		{name='S', code=self:template[[	value.v<?=vec3?> = <?=vec3?>_cross(calc_E(*U), calc_H(*U));]], type=env.vec3, units='kg/s^3'},
		{
			name='energy', code=self:template[[
	<?=susc_t?> _1_eps = <?=susc_t?>_mul(U->sqrt_1_eps, U->sqrt_1_eps);
	<?=susc_t?> _1_mu = <?=susc_t?>_mul(U->sqrt_1_mu, U->sqrt_1_mu);
	value.vreal = <?=real_mul?>(<?=add?>(
		<?=scalar?>_<?=susc_t?>_mul(eqn_coordLenSq(U->D, x), _1_eps),
		<?=scalar?>_<?=susc_t?>_mul(eqn_coordLenSq(calc_H(*U), x), _1_mu)
	), .5);
]],
			type = scalar,
			units = 'kg/(m*s^2)',
		},
	}

	vars:insert(self:createDivDisplayVar{field='D', scalar=env.scalar, units='C/m^3'} or nil)
	vars:insert(self:createCurlDisplayVar{field='D', scalar=env.scalar, units='C/m^3'} or nil)
	vars:insert(self:createDivDisplayVar{field='B', scalar=env.scalar, units='kg/(C*m*s)'} or nil)
	vars:insert(self:createCurlDisplayVar{field='B', scalar=env.scalar, units='kg/(C*m*s)'} or nil)

	return vars
end

function Maxwell:eigenWaveCodePrefix(n, eig, x, waveIndex)
--[=[
	return self:template([[
	<?=scalar?> v_p_abs = <?=mul?>(<?=eig?>.sqrt_1_eps, <?=eig?>.sqrt_1_mu);
]], {
		eig = '('..eig..')',
	})
--]=]
-- [=[
	local env = self:getEnv()
	local code = self:template(
		[[<?=mul?>(<?=eig?>.sqrt_1_eps, <?=eig?>.sqrt_1_mu)]],
		{
			eig = '('..eig..')',
		}
	)
	if self.scalar == 'cplx' then
		code = env.abs..'('..code..')'
	end
	return 'real v_p_abs = '..code..';'
--]=]
end

function Maxwell:eigenWaveCode(n, eig, x, waveIndex)
	waveIndex = math.floor(waveIndex / self.numRealsInScalar)
	return ({
		'-v_p_abs',
		'-v_p_abs',
		'0',
		'0',
		'v_p_abs',
		'v_p_abs',
	})[waveIndex+1] or error('got a bad waveIndex: '..waveIndex)
end

function Maxwell:eigenMaxWaveCode(n, eig, x)
	return 'v_p_abs'
end
function Maxwell:eigenMinWaveCode(n, eig, x)
	return '-'..self:eigenMaxWaveCode(n, eig, x)
end

function Maxwell:consWaveCodePrefix(n, U, x, waveIndex)
	local env = self:getEnv()
	local code = self:template(
		[[<?=mul?>(<?=U?>.sqrt_1_eps, <?=U?>.sqrt_1_mu)]],
		{
			U = '('..U..')',
		}
	)
	if self.scalar == 'cplx' then
		code = env.abs..'('..code..')'
	end
	return 'real v_p_abs = '..code..';'
end
Maxwell.consWaveCode = Maxwell.eigenWaveCode

function Maxwell:consMaxWaveCode(n, U, x)
	return 'v_p_abs'
end
function Maxwell:consMinWaveCode(n, U, x)
	return '-'..self:consMaxWaveCode(n, U, x)
end

return Maxwell

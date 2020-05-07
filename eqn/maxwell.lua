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
local range = require 'ext.range'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local template = require 'template'
local common = require 'common'
local xNames = common.xNames

local Maxwell = class(Equation)
Maxwell.name = 'Maxwell'

Maxwell.hasEigenCode = true
Maxwell.hasFluxFromConsCode = true
Maxwell.roeUseFluxFromCons = true
Maxwell.useSourceTerm = true

-- don't incorporate the Conn^k_ij E_k terms into the flux
Maxwell.weightFluxByGridVolume = false

Maxwell.initStates = require 'init.euler'

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
		{name='sigma', type=self.scalar, units='(C^2*s)/(kg*m^3)'},
		{name='sqrt_1_eps', type=self.susc_t, units='(kg*m^3)^.5/(C*s)'},
		{name='sqrt_1_mu', type=self.susc_t, units='C/(kg*m)^.5'},
	}
	
	self.eigenVars = table{
		{name='sqrt_1_eps', type=self.susc_t, units='(kg*m^3)^.5/(C*s)'},
		{name='sqrt_1_mu', type=self.susc_t, units='C/(kg*m)^.5'},
	}

	Maxwell.super.init(self, args)

	self.postComputeFluxCode = template([[
<? local vec3 = eqn.vec3 ?>
		//TODO shouldn't I be transforming both the left and right fluxes by the metrics at their respective coordinates?
		//flux is computed raised via Levi-Civita upper
		//so here we lower it
		real _1_sqrt_det_g = 1. / coord_sqrt_det_g(x);
		flux.D = <?=vec3?>_real_mul(eqn_coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = <?=vec3?>_real_mul(eqn_coord_lower(flux.B, x), _1_sqrt_det_g);
]], {eqn=self})


	local NoDiv = require 'op.nodiv'{
		poissonSolver = require 'op.poisson_jacobi',
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

function Maxwell:getCommonFuncCode()
	return template([[
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

]], self:getTemplateEnv())
end

Maxwell.initStateCode = [[

<? 
local cons_t = eqn.cons_t
local susc_t = eqn.susc_t
local scalar = eqn.scalar
local vec3 = eqn.vec3
local zero = scalar..'_zero'
?>

kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
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
	U->sigma = conductivity;
	U->rhoCharge = <?=zero?>;
	U->sqrt_1_eps = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permittivity));
	U->sqrt_1_mu = <?=susc_t?>_sqrt(<?=susc_t?>_inv(permeability));
}
]]

function Maxwell:getSolverCode()
	return template(file['eqn/maxwell.cl'], self:getTemplateEnv())
end

function Maxwell:getTemplateEnv()
	local scalar = self.scalar
	local env = {}
	env.eqn = self
	env.solver = self.solver
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

-- k is 0,1,2
local function curl(eqn,k,result,field,env)
	local i = (k+1)%3
	local j = (i+1)%3
	return {['curl '..field..' '..xNames[k+1]] = template([[
	if (OOB(1,1)) {
		<?=result?> = <?=zero?>;
	} else {

<? if i+1 <= solver.dim then ?>
		global const <?=eqn.cons_t?>* Uim = U - solver->stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Uip = U + solver->stepsize.s<?=i?>;
		<?=scalar?> vim_j = Uim-><?=field?>.s<?=j?>;
		<?=scalar?> vip_j = Uip-><?=field?>.s<?=j?>;
<? else ?>
		<?=scalar?> vim_j = <?=zero?>;
		<?=scalar?> vip_j = <?=zero?>;
<? end?>

<? if j+1 <= solver.dim then ?>
		global const <?=eqn.cons_t?>* Ujm = U - solver->stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + solver->stepsize.s<?=j?>;
		<?=scalar?> vjm_i = Ujm-><?=field?>.s<?=i?>;
		<?=scalar?> vjp_i = Ujp-><?=field?>.s<?=i?>;
<? else ?>
		<?=scalar?> vjm_i = <?=zero?>;
		<?=scalar?> vjp_i = <?=zero?>;
<? end ?>

		<?=result?> = <?=sub?>(
			<?=real_mul?>(<?=sub?>(vjp_i, vjm_i), 1. / (2. * solver->grid_dx.s<?=i?>)),
			<?=real_mul?>(<?=sub?>(vip_j, vim_j), 1. / (2. * solver->grid_dx.s<?=j?>))
		);
	}
]], table(env, {
		i = i,
		j = j,
		result = result,
		field = field,
	}))}
end

Maxwell.predefinedDisplayVars = {
	'U D x',
	'U D y',
	'U D z',
	'U D mag',
	'U div D',
	'U phi',
	'U B x',
	'U B y',
	'U B z',
	'U B mag',
	'U div B',
	'U psi',
}

function Maxwell:getDisplayVars()
	local env = self:getTemplateEnv()
	
	local vars = Maxwell.super.getDisplayVars(self)
	vars:append{ 
		{name='E', code=template([[	value.v<?=vec3?> = calc_E(*U);]], env), type=env.vec3, units='(kg*m)/(C*s)'},
		{name='H', code=template([[	value.v<?=vec3?> = calc_H(*U);]], env), type=env.vec3, units='C/(m*s)'},
		{name='S', code=template([[	value.v<?=vec3?> = <?=vec3?>_cross(calc_E(*U), calc_H(*U));]], env), type=env.vec3, units='kg/s^3'},
		{
			name='energy', code=template([[
	<?=susc_t?> _1_eps = <?=susc_t?>_mul(U->sqrt_1_eps, U->sqrt_1_eps);
	<?=susc_t?> _1_mu = <?=susc_t?>_mul(U->sqrt_1_mu, U->sqrt_1_mu);
	value.vreal = <?=real_mul?>(<?=add?>(
		<?=scalar?>_<?=susc_t?>_mul(eqn_coordLenSq(U->D, x), _1_eps),
		<?=scalar?>_<?=susc_t?>_mul(eqn_coordLenSq(calc_H(*U), x), _1_mu)
	), .5);
]], env), 
			type = scalar,
			units = 'kg/(m*s^2)',
		},
	}:append(table{'D', 'B'}:map(function(field,i)
		local field = assert( ({D='D', B='B'})[field] )
		return {
			name='div '..field, 
			code=template([[
	<?=scalar?> v = <?=zero?>;
<? for j=0,solver.dim-1 do ?>
	v = <?=add?>(v, <?=real_mul?>(
		<?=sub?>(
			U[solver->stepsize.s<?=j?>].<?=field?>.s<?=j?>,
			U[-solver->stepsize.s<?=j?>].<?=field?>.s<?=j?>
		), .5 / solver->grid_dx.s<?=j?>));
<? end ?>
	value.v<?=scalar?> = v;
]], table(env, {field=field})),
			type = scalar, 
			units = ({
				D = 'C/m^3',
				B = 'kg/(C*m*s)',
			})[field],
		}
	end))

	for _,field in ipairs{'D', 'B'} do
		local v = range(0,2):map(function(i)
			return curl(self,i,'value.v'..env.vec3..'.s'..i,field, env)
		end)
		vars:insert{
			name='curl '..field, 
			code=template([[
	<? for i=0,2 do ?>{
		<?=select(2,next(v[i+1]))?>
	}<? end ?>
]], {v=v}),
			type = env.vec3,
			units = ({
				D = 'C/m^3',
				B = 'kg/(C*m*s)',
			})[field],
		}
	end

	return vars
end

function Maxwell:eigenWaveCodePrefix(n, eig, x, waveIndex)
--[=[
	return template([[
	<?=scalar?> v_p_abs = <?=mul?>(<?=eig?>.sqrt_1_eps, <?=eig?>.sqrt_1_mu);
]], table(self:getTemplateEnv(), {
		eqn = self,
		eig = '('..eig..')',
	}))
--]=]
-- [=[
	local env = self:getTemplateEnv()
	local code = template(
		[[<?=mul?>(<?=eig?>.sqrt_1_eps, <?=eig?>.sqrt_1_mu)]],
		table(env, {
			eqn = self,
			eig = '('..eig..')',
		})
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
	local env = self:getTemplateEnv()
	local code = template(
		[[<?=mul?>(<?=U?>.sqrt_1_eps, <?=U?>.sqrt_1_mu)]],
		table(env, {
			eqn = self,
			U = '('..U..')',
		})
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

--[[ 2000 Munz
But then I changed the eigendecomposition to match my Trangenstein one with epsilons and mus
And then I changed epsilon E -> D, so it's based on B and D 

I'm still not completely sold on this technique.
Relativistically we are doing:
~F^uv = F^uv + eta^uv phi chi, ~F^uv_,v = mu j^u
~*F^uv = *F^uv + eta^uv psi gamma, ~*F^uv_,v = 0
except with the special twist of changing the coefficients next to the eta^tt term from chi or gamma to 1/chi or 1/gamma
I'm still suspicious on this special exception.
Not doing so means our wavespeeds are all the speed of light, and we can't control the propagation of GLM variables.

--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local clnumber = require 'cl.obj.number'
local template = require 'template'

local GLM_Maxwell = class(Equation)
GLM_Maxwell.name = 'GLM_Maxwell'

GLM_Maxwell.mirrorVars = {{'D.x', 'B.x'}, {'D.y', 'B.y'}, {'D.z', 'B.z'}}

GLM_Maxwell.hasEigenCode = true
GLM_Maxwell.hasFluxFromConsCode = true
GLM_Maxwell.useSourceTerm = true
GLM_Maxwell.roeUseFluxFromCons = true

GLM_Maxwell.initStates = require 'init.euler'


function GLM_Maxwell:init(args)

	self.scalar = 'real'
	--self.scalar = 'cplx'
	self.vec3 = self.scalar..'3'

	-- TODO tensor susceptibilty support ... but that affects the eigendecomposition ...
	self.susc_t = self.scalar

	self.numRealsInScalar = ffi.sizeof(self.scalar) / ffi.sizeof'real'

	self.numIntStates = 8 * self.numRealsInScalar
	self.numWaves = 8 * self.numRealsInScalar


	self.consVars = {
		{D = self.vec3},
		{B = self.vec3},
		{phi = self.scalar},
		{psi = self.scalar},
		{rhoCharge = self.scalar},
		{sigma = self.scalar},
		{sqrt_1_eps = self.susc_t},
		{sqrt_1_mu = self.susc_t},
	}

	self.eigenVars = table{
		{sqrt_1_eps = self.scalar},
		{sqrt_1_mu = self.scalar},
	}


	GLM_Maxwell.super.init(self, args)


	self.postComputeFluxCode = template([[
<? local vec3 = eqn.vec3 ?>
		//TODO shouldn't I be transforming both the left and right fluxes by the metrics at their respective coordinates?
		//flux is computed raised via Levi-Civita upper
		//so here we lower it
		real _1_sqrt_det_g = 1. / sqrt_det_g_grid(x);
		flux.D = <?=vec3?>_real_mul(eqn_coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = <?=vec3?>_real_mul(eqn_coord_lower(flux.B, x), _1_sqrt_det_g);
	]], {eqn=self})
end

function GLM_Maxwell:getCommonFuncCode()
	return template([[
<? -- in common with Maxwell ?>
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
	return <?=vec3?>_<?=scalar?>_mul(U.D, <?=mul?>(U.sqrt_1_eps, U.sqrt_1_eps));
}
<?=vec3?> calc_H(<?=eqn.cons_t?> U) { 
	return <?=vec3?>_<?=scalar?>_mul(U.B, <?=mul?>(U.sqrt_1_mu, U.sqrt_1_mu));
}

real ESq(<?=eqn.cons_t?> U, real3 x) { return eqn_coordLenSq(calc_E(U), x); }
real BSq(<?=eqn.cons_t?> U, real3 x) { return eqn_coordLenSq(U.B, x); }
]], self:getTemplateEnv())
end

GLM_Maxwell.initStateCode = [[

<? 
local cons_t = eqn.cons_t
local susc_t = eqn.susc_t
local scalar = eqn.scalar
local vec3 = eqn.vec3
local zero = scalar..'_zero'
local inv = scalar..'_inv'
local fromreal = scalar..'_from_real'
local sqrt = scalar..'_sqrt'
?>

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
	<?=vec3?> D = <?=vec3?>_zero;
	<?=vec3?> B = <?=vec3?>_zero;
	<?=scalar?> conductivity = <?=fromreal?>(1.);
	
	<?=scalar?> permittivity = <?=fromreal?>(1.);
	<?=scalar?> permeability = <?=fromreal?>(1.);
	
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
	U->sqrt_1_eps = <?=sqrt?>(<?=inv?>(permittivity));
	U->sqrt_1_mu = <?=sqrt?>(<?=inv?>(permeability));
}
]]

function GLM_Maxwell:getSolverCode()
	return template(file['eqn/glm-maxwell.cl'], self:getTemplateEnv())
end

function GLM_Maxwell:getTemplateEnv()
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
	local xs = {'x','y','z'}
	local i = (k+1)%3
	local j = (i+1)%3
	return {['curl '..field..' '..xs[k+1]] = template([[
	if (OOB(1,1)) {
		<?=result?> = <?=zero?>;
	} else {

<? if i+1 <= solver.dim then ?>
		global const <?=eqn.cons_t?>* Uim = U - stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Uip = U + stepsize.s<?=i?>;
		<?=scalar?> vim_j = Uim-><?=field?>.s<?=j?>;
		<?=scalar?> vip_j = Uip-><?=field?>.s<?=j?>;
<? else ?>
		<?=scalar?> vim_j = <?=zero?>;
		<?=scalar?> vip_j = <?=zero?>;
<? end?>

<? if j+1 <= solver.dim then ?>
		global const <?=eqn.cons_t?>* Ujm = U - stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + stepsize.s<?=j?>;
		<?=scalar?> vjm_i = Ujm-><?=field?>.s<?=i?>;
		<?=scalar?> vjp_i = Ujp-><?=field?>.s<?=i?>;
<? else ?>
		<?=scalar?> vjm_i = <?=zero?>;
		<?=scalar?> vjp_i = <?=zero?>;
<? end ?>

		<?=result?> = <?=sub?>(
			<?=real_mul?>(<?=sub?>(vjp_i, vjm_i), 1. / (2. * grid_dx<?=i?>)),
			<?=real_mul?>(<?=sub?>(vip_j, vim_j), 1. / (2. * grid_dx<?=j?>))
		);
	}
]], table(env, {
		i = i,
		j = j,
		result = result,
		field = field,
	}))}
end

function GLM_Maxwell:getDisplayVars()
	local env = self:getTemplateEnv()
	
	local vars = GLM_Maxwell.super.getDisplayVars(self)
	vars:append{ 
		{S = template([[
	*value_<?=vec3?> = <?=vec3?>_cross(calc_E(*U), calc_H(*U));
]], env), type='real3'},
		{energy = template([[
	*value = (eqn_coordLenSq(U->D, x) + eqn_coordLenSq(calc_H(*U), x)) * .5;
]], env)},
	}:append(table{'D','B'}:map(function(var,i)
		local field = assert( ({D='D', B='B'})[var] )
		return {['div '..var] = template([[
	<?=scalar?> v = <?=zero?>;
<? for j=0,solver.dim-1 do ?>
	v = <?=add?>(v, <?=real_mul?>(
		<?=sub?>(
			U[stepsize.s<?=j?>].<?=field?>.s<?=j?>,
			U[-stepsize.s<?=j?>].<?=field?>.s<?=j?>
	), 1. / grid_dx<?=j?>));
<? end ?>	
	v = <?=real_mul?>(v, .5);
	*value_<?=scalar?> = v;
]], table(env, {solver=self.solver, field=field}))}
	end))

	for _,field in ipairs{'D', 'B'} do
		local v = range(0,2):map(function(i) 
			return curl(self,i,'value_'..env.vec3..'->s'..i,field, env) 
		end)
		vars:insert{['curl '..field]= template([[
	<? for i=0,2 do ?>{
		<?=select(2,next(v[i+1]))?>
	}<? end ?>
]], {v=v}), type='real3'}
	end

	return vars
end

function GLM_Maxwell:eigenWaveCodePrefix(side, eig, x, waveIndex)
--[=[	
	return template([[
	<?=scalar?> v_p = <?=mul?>(<?=eig?>.sqrt_1_eps, <?=eig?>.sqrt_1_mu);
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


-- to use this, I really need cplx multiplications everywhere it is used
-- which is in the roe solver and hll solver
function GLM_Maxwell:eigenWaveCode(side, eig, x, waveIndex)
	waveIndex = math.floor(waveIndex / self.numRealsInScalar)
	return template(({
		'-v_p_abs * divPhiWavespeed',
		'-v_p_abs * divPsiWavespeed',
		'-v_p_abs',
		'-v_p_abs',
		'v_p_abs',
		'v_p_abs',
		'v_p_abs * divPhiWavespeed',
		'v_p_abs * divPsiWavespeed',
	})[waveIndex+1] or error('got a bad waveIndex: '..waveIndex), 
		self:getTemplateEnv())
end

function GLM_Maxwell:eigenMaxWaveCode(side, eig, x)
	return 'max(max(divPsiWavespeed, divPhiWavespeed), 1.) * v_p_abs;'
end
function GLM_Maxwell:eigenMinWaveCode(side, eig, x)
	return '-'..self:eigenMaxWaveCode(side, eig, x)
end


function GLM_Maxwell:consWaveCodePrefix(side, U, x, waveIndex) 
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
GLM_Maxwell.consWaveCode = GLM_Maxwell.eigenWaveCode

function GLM_Maxwell:consMaxWaveCode(side, U, x)
	return 'max(max(divPsiWavespeed, divPhiWavespeed), 1.) * v_p_abs;'
end
function GLM_Maxwell:consMinWaveCode(side, U, x)
	return '-'..self:consMaxWaveCode(side, U, x)
end


return GLM_Maxwell

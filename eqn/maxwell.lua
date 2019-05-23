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
local file = require 'ext.file'
local range = require 'ext.range'
local Equation = require 'eqn.eqn'
local template = require 'template'
local common = require 'common'()
local xNames = common.xNames

local Maxwell = class(Equation)

-- don't incorporate the Conn^k_ij E_k terms into the flux
Maxwell.weightFluxByGridVolume = false

Maxwell.name = 'Maxwell'
Maxwell.numIntStates = 6

-- I'm working on making complex numbers exchangeable
Maxwell.scalar = 'real'
--Maxwell.scalar = 'cplx'

-- 1 for 'real', 2 for 'cplx'
Maxwell.numRealsInScalar = ffi.sizeof(Maxwell.scalar) / ffi.sizeof'real'

Maxwell.numWaves = 6 * Maxwell.numRealsInScalar

Maxwell.vec3 = Maxwell.scalar..'3'
Maxwell.mat3x3 = Maxwell.scalar..'3x3'

Maxwell.consVars = table{
	{D = Maxwell.vec3},				-- C/m^2
	{B = Maxwell.vec3},				-- kg/(C s)
	
	{DPot = Maxwell.scalar},		-- C/m^2
	{BPot = Maxwell.scalar},		-- kg/(C s)
	
	{rhoCharge = Maxwell.scalar},	-- C/m^3
	{sigma = Maxwell.scalar},		-- (C^2 s)/(kg m^3)
}

--[[
TODO make these complex
but that means making E and B complex 
and that means complex math, and *drumroll* complex code generation of the coordLenSq functions
and this would be easier if OpenCL supported the 'complex' keyword

another todo - max this a tensor
some common susceptibility tensors are symmetric?
I thought I caught somewhere that they are often projection matrices...
https://physics.stackexchange.com/questions/351012/how-can-i-deduce-the-magnetic-susceptibility-tensor-of-a-biaxial-liquid-crystal
https://physics.stackexchange.com/questions/148634/equivalent-tensor-order-parameters-of-nematic-liquid-crystals?rq=1
https://en.wikipedia.org/wiki/Permittivity#Tensorial_permittivity
https://en.wikipedia.org/wiki/Electro-gyration
--]]

Maxwell.susc_t = Maxwell.scalar

--[[ 
there's a catch, if susceptibility is a tensor, 
then it gets applied to the flux matrix, 
and therefore it needs to be considered in the flux decomposition

Therefore I don't just need the susceptibility tensor,
I also need the eigen-decomposition of the susceptibility times the Levi-Civita tensor

So I need to store the U, S, and V matrices of this...
So we need three matrices and not just one ...
--]]
--Maxwell.susc_t = Maxwell.mat3x3

Maxwell.consVars:append{
	{sqrt_1_eps = Maxwell.susc_t},		-- sqrt( (kg m^3)/(C^2 s^2) )
	{sqrt_1_mu = Maxwell.susc_t},		-- sqrt( C^2/(kg m) )
}

Maxwell.mirrorVars = {{'D.x', 'B.x'}, {'D.y', 'B.y'}, {'D.z', 'B.z'}}

Maxwell.hasEigenCode = true
Maxwell.hasFluxFromConsCode = true
Maxwell.useSourceTerm = true
Maxwell.roeUseFluxFromCons = true

Maxwell.initStates = require 'init.euler'

Maxwell.postComputeFluxCode = template([[
<? local vec3 = eqn.vec3 ?>
		//TODO shouldn't I be transforming both the left and right fluxes by the metrics at their respective coordinates?
		//flux is computed raised via Levi-Civita upper
		//so here we lower it
		real _1_sqrt_det_g = 1. / sqrt_det_g_grid(x);
		flux.D = <?=vec3?>_real_mul(eqn_coord_lower(flux.D, x), _1_sqrt_det_g);
		flux.B = <?=vec3?>_real_mul(eqn_coord_lower(flux.B, x), _1_sqrt_det_g);
]], {eqn=Maxwell})

function Maxwell:init(args)
	Maxwell.super.init(self, args)

	local NoDiv = class(require 'op.nodiv')
	NoDiv.scalar = self.scalar
	self.solver.ops:insert(NoDiv{
		solver = self.solver,
	})
	-- should I be fixing div E = rhoCharge, 
	-- or should I get rid of the rhoCharge field and the div E constraint?
	self.solver.ops:insert(NoDiv{
		solver = self.solver,
		potentialField = 'DPot',
		chargeField = 'rhoCharge',
	})
end

function Maxwell:getCommonFuncCode()
	return template([[
//hmm, for E and B, even if the coord is 2D, we need all 3D components ...
//this means we need coordLen functions with guaranteed dimensions, including tangent spaces

<? if eqn.scalar == 'real' then ?>

real eqn_coordLenSq(real3 v, real3 x) {
	return coordLenSq(v, x);
}

#define eqn_cartesianToCoord cartesianToCoord
#define eqn_coord_lower coord_lower

<? elseif eqn.scalar == 'cplx' then ?>

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

<? end -- eqn.scalar ?>

<?=eqn.vec3?> calc_E(<?=eqn.cons_t?> U) { 
/*
scalar type / susceptibility type:
real/scalar => real_real3_mul
cplx/scalar => cplx_cplx3_mul
real/tensor => real3x3_real3_mul
cplx/tensor => cplx3x3_cplx3_mul
*/
	return <?=eqn.susc_t?>_<?=eqn.vec3?>_mul(
		<?=eqn.susc_t?>_mul(U.sqrt_1_eps, U.sqrt_1_eps),
		U.D);
}

<?=eqn.vec3?> calc_H(<?=eqn.cons_t?> U) { 
	return <?=eqn.susc_t?>_<?=eqn.vec3?>_mul(
		<?=eqn.susc_t?>_mul(U.sqrt_1_mu, U.sqrt_1_mu),
		U.B);
}
/*
|E| = E_i *E^i = (re E^i + i im E^i)  (re E^j - i im E^j) gamma_ij
= ((re E^i re E^j + im E^i im E^j) + i (im E^i re E^j - re E^i im E^j)) gamma_ij
= (re E^i re E^j + im E^i im E^j) gamma_ij
re |E| = coordLenSq(re_E, re_E) + coordLenSq(im_E, im_E)
im |E| = 0
*/
real ESq(<?=eqn.cons_t?> U, real3 x) {
	return eqn_coordLenSq(calc_E(U), x);
}

real BSq(<?=eqn.cons_t?> U, real3 x) {
	return eqn_coordLenSq(U.B, x);
}

]], {
		eqn = self,
	})
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
	global <?=solver.solver_t?>* solver,
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
	
	U->D = <?=vec3?>_real_mul(eqn_cartesianToCoord(D, x), unit_C_per_m2);
	U->B = <?=vec3?>_real_mul(eqn_cartesianToCoord(B, x), unit_kg_per_C_s);
	U->BPot = <?=zero?>;
	U->sigma = <?=susc_t?>_real_mul(conductivity, unit_C2_s_per_kg_m3);
	U->sqrt_1_eps = <?=susc_t?>_sqrt(<?=susc_t?>_inv(<?=susc_t?>_real_mul(permittivity, unit_C2_s2_per_kg_m3)));
	U->sqrt_1_mu = <?=susc_t?>_sqrt(<?=susc_t?>_inv(<?=susc_t?>_real_mul(permeability, unit_kg_m_per_C2)));
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
	env.real_mul = scalar..'_real_mul'
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

--[[
for E = [0, sin(x-t), 0]
dEy/dx = cos(x-t)
so curl(E).z = -cos(x-t)
--]]
function Maxwell:getDisplayVars()
	local env = self:getTemplateEnv()

	local vars = Maxwell.super.getDisplayVars(self)
	vars:append{ 
		{E = template([[	*value_<?=vec3?> = calc_E(*U);]], env), type=vec3},
		{S = template([[
	*value_<?=vec3?> = <?=vec3?>_cross(calc_E(*U), calc_H(*U));
]], env), type=vec3},
		{energy = template([[
	*value = (eqn_coordLenSq(U->D, x) + eqn_coordLenSq(calc_H(*U), x)) * .5;
]], env), type=scalar},
	}:append(table{'D','B'}:map(function(var,i)
		local field = assert( ({D='D', B='B'})[var] )
		return {['div '..var] = template([[
	<?=scalar?> v = <?=zero?>;
<? for j=0,solver.dim-1 do ?>
	v = <?=add?>(v, <?=real_mul?>(
		<?=sub?>(
			U[solver->stepsize.s<?=j?>].<?=field?>.s<?=j?>,
			U[-solver->stepsize.s<?=j?>].<?=field?>.s<?=j?>
		), 1. / solver->grid_dx.s<?=j?>));
<? end ?>

	v = <?=real_mul?>(v, .5);
	*value_<?=scalar?> = v;
]], table(env, {field=field}))}
	end))

	for _,field in ipairs{'D', 'B'} do
		local v = range(0,2):map(function(i)
			return curl(self,i,'value_'..env.vec3..'->s'..i,field,env)
		end)
		vars:insert{['curl '..field]= template([[
	<? for i=0,2 do ?>{
		<?=select(2,next(v[i+1]))?>
	}<? end ?>
]], {v=v}), type=vec3}
	end

	return vars 
end

Maxwell.eigenVars = table{
	{sqrt_1_eps = Maxwell.susc_t},
	{sqrt_1_mu = Maxwell.susc_t},
}

function Maxwell:eigenWaveCodePrefix(side, eig, x, waveIndex)
	return template([[
	<?=eqn.susc_t?> eig_lambda = <?=eqn.susc_t?>_mul(<?=eig?>.sqrt_1_eps, <?=eig?>.sqrt_1_mu);
]], {
		eqn = self,
		eig = '('..eig..')',
	})
end

function Maxwell:eigenWaveCode(side, eig, x, waveIndex)
	if self.susc_t == 'real' then
		if waveIndex == 0 or waveIndex == 1 then
			return '-eig_lambda'
		elseif waveIndex == 2 or waveIndex == 3 then
			return '0'
		elseif waveIndex == 4 or waveIndex == 5 then
			return 'eig_lambda'
		else
			error'got a bad waveIndex'
		end
	elseif self.susc_t == 'cplx' then
		if waveIndex == 0 or waveIndex == 2 then
			return '-eig_lambda.re'
		elseif waveIndex == 1 or waveIndex == 3 then
			return '-eig_lambda.im'
		elseif waveIndex >= 4 or waveIndex <= 7 then
			return '0'
		elseif waveIndex == 8 or waveIndex == 10 then
			return 'eig_lambda.re'
		elseif waveIndex == 9 or waveIndex == 11 then
			return 'eig_lambda.im'
		else
			error'got a bad waveIndex'
		end
	end
end

function Maxwell:consWaveCodePrefix(side, U, x, waveIndex)
	return template([[
<? local susc_t = eqn.susc_t ?>
	<?=susc_t?> eig_lambda = <?=susc_t?>_mul(<?=U?>.sqrt_1_eps, <?=U?>.sqrt_1_mu);
]], {
		eqn = self,
		U = '('..U..')',
	})
end
Maxwell.consWaveCode = Maxwell.eigenWaveCode

return Maxwell

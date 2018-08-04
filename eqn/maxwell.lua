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
	{D = Maxwell.vec3},
	{B = Maxwell.vec3},
	
	{DPot = Maxwell.scalar},
	{BPot = Maxwell.scalar},
	
	{rhoCharge = Maxwell.scalar},
	{sigma = Maxwell.scalar},
}
	
--[[
TODO make these complex
but that means making E and B complex 
and that means complex math, and *drumroll* complex code generation of the coordLenSq functions
and this would be easier if OpenCL supported the 'complex' keyword

another todo - max this a tensor
some common susceptibility tensors are symmetric?  I thought I caught somewhere that they are often projection matrices...
--]]

Maxwell.susceptibilityTensor = false

if Maxwell.susceptibilityTensor then
	Maxwell.consVars:append{
		{_1_eps = Maxwell.mat3x3},
		{_1_mu = Maxwell.mat3x3},
	}
else
	Maxwell.consVars:append{
		{_1_eps = Maxwell.scalar},
		{_1_mu = Maxwell.scalar},
	}
end

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

cplx3 eqn_cartesianToCoord(real3 v, real3 x) {
	return cplx3_from_real3(
		cartesianToCoord(v, x),
		cartesianToCoord(_real3(0,0,0), x));
}

cplx3 eqn_coord_lower(cplx3 v, real3 x) {
	return cplx3_from_real3(
		coord_lower(cplx3_re(v), x),
		coord_lower(cplx3_im(v), x));
}

<? end -- eqn.scalar ?>

<?=eqn.vec3?> calc_E(<?=eqn.cons_t?> U) {
	return <?=eqn.vec3?>_scale(U.D, U._1_eps);
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

<? local fromreal = eqn.scalar..'_from_real' ?>

kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
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
	real3 E = _real3(0,0,0);
	real3 B = _real3(0,0,0);
	real conductivity = 1.;
	
	real permittivity = 1.;
	real permeability = 1.;
	
	//throw-away
	real rho = 0;
	real3 v = _real3(0,0,0);
	real P = 0;
	real ePot = 0;
	
	<?=code?>
	
	U->D = eqn_cartesianToCoord(real3_scale(E, permittivity), x);
	U->B = eqn_cartesianToCoord(B, x);
	U->BPot = <?=fromreal?>(0);
	U->sigma = <?=fromreal?>(conductivity);
	U->_1_eps = <?=fromreal?>(1. / permittivity);
	U->_1_mu = <?=fromreal?>(1. / permeability);
}
]]

function Maxwell:getSolverCode()
	return template(file['eqn/maxwell.cl'], {eqn=self})
end

-- k is 0,1,2
local function curl(eqn,k,result,field)
	local i = (k+1)%3
	local j = (i+1)%3
	return {['curl '..field..' '..xNames[k+1]] = template([[
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

--[[
for E = [0, sin(x-t), 0]
dEy/dx = cos(x-t)
so curl(E).z = -cos(x-t)
--]]
function Maxwell:getDisplayVars()
	if self.scalar == 'real' then
		local vars = Maxwell.super.getDisplayVars(self)
		vars:append{ 
			{E = '*valuevec = real3_scale(U->D, U->_1_eps);', type='real3'},
			{S = '*valuevec = real3_scale(real3_cross(U->D, U->B), U->_1_eps);', type='real3'},
			{energy = [[
	*value = .5 * (coordLenSq(U->D, x) + coordLenSq(U->B, x) * U->_1_mu * U->_1_mu);
]]},
		}:append(table{'E','B'}:map(function(var,i)
			local field = assert( ({E='D', B='B'})[var] )
			return {['div '..var] = template([[
	*value = .5 * (0.
<?
for j=0,solver.dim-1 do
?>		+ (U[stepsize.s<?=j?>].<?=field?>.s<?=j?> 
			- U[-stepsize.s<?=j?>].<?=field?>.s<?=j?>
		) / grid_dx<?=j?>
<?
end
?>	)<? 
if field == 'D' then 
?> * U->_1_eps<?
end
?>;
]], {solver=self.solver, field=field})}
		end))

		for _,field in ipairs{'D', 'B'} do
			local v = range(0,2):map(function(i) 
				return curl(self,i,'valuevec->s'..i,field) 
			end)
			vars:insert{['curl '..field]= template([[
		<? for i=0,2 do ?>{
			<?=select(2,next(v[i+1]))?>
		}<? end ?>
	]], {v=v}), type='real3'}
		end

		return vars 
	elseif self.scalar == 'cplx' then
		local vars = {
			{['D re'] = '*valuevec = cplx3_re(U->D);', type='real3'},
			{['D im'] = '*valuevec = cplx3_im(U->D);', type='real3'},
			{['B re'] = '*valuevec = cplx3_re(U->B);', type='real3'},
			{['B im'] = '*valuevec = cplx3_im(U->B);', type='real3'},
		}
		return vars
	end
end

Maxwell.eigenVars = table{
	{sqrt_1_eps = Maxwell.scalar},
	{sqrt_1_mu = Maxwell.scalar},
}

function Maxwell:eigenWaveCodePrefix(side, eig, x, waveIndex)
	return template([[
	<?=eqn.scalar?> eig_lambda = <?=eqn.scalar?>_mul(<?=eig?>.sqrt_1_eps, <?=eig?>.sqrt_1_mu);
]], {
		eqn = self,
		eig = '('..eig..')',
	})
end

function Maxwell:eigenWaveCode(side, eig, x, waveIndex)
	if self.scalar == 'real' then
		if waveIndex == 0 or waveIndex == 1 then
			return '-eig_lambda'
		elseif waveIndex == 2 or waveIndex == 3 then
			return '0'
		elseif waveIndex == 4 or waveIndex == 5 then
			return 'eig_lambda'
		else
			error'got a bad waveIndex'
		end
	elseif self.scalar == 'cplx' then
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
<? local scalar = eqn.scalar ?>
	<?=scalar?> eig_lambda = <?=scalar?>_sqrt(<?=scalar?>_mul(<?=U?>._1_eps, <?=U?>._1_mu));
]], {
		eqn = self,
		U = '('..U..')',
	})
end
Maxwell.consWaveCode = Maxwell.eigenWaveCode

return Maxwell

local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'


local Euler = class(Equation)
Euler.name = 'Euler'

-- ePot is the 6th param
-- which means it's now in the derivBuf, but it is always zero
-- so TODO a new variable for deriv size vs cons_t size?
Euler.numStates = 6	

Euler.numWaves = 5
Euler.numIntStates = 5	-- don't bother integrate ePot

Euler.mirrorVars = {{'m.x'}, {'m.y'}, {'m.z'}} 
Euler.hasCalcDTCode = true
Euler.hasEigenCode = true
Euler.hasFluxFromConsCode = true
Euler.roeUseFluxFromCons = true

-- the only source term that the Euler equations has is the connection coefficients of the velocity vector
-- maybe later I will automatically flag what elements are vectors
-- and automatically add connection coefficients
--
-- I'm also using the source term for the viscosity ... if you choose to use FANS
Euler.useSourceTerm = true

Euler.initStates = require 'init.euler'

function Euler:init(args)
	Euler.super.init(self, args)

	local SelfGrav = require 'op.selfgrav'

if require 'solver.meshsolver'.is(self.solver) then
	print("not using selfgrav with mesh solvers yet")
else
		self.gravOp = SelfGrav{solver = self.solver}
		self.solver.ops:insert(self.gravOp)
end

end

Euler.primVars = table{
	{rho = 'real'},
	{v = 'real3'},
	{P = 'real'},
	{ePot = 'real'},
}

Euler.consVars = table{
	{rho = 'real'},
	{m = 'real3'},
	{ETotal = 'real'},
	{ePot = 'real'},
}

function Euler:createInitState()
	Euler.super.createInitState(self)
	self:addGuiVars{
		{name='heatCapacityRatio', value=7/5},
	}
end

function Euler:getCommonFuncCode()
	return template([[
inline real calc_H(real P) { return P * (heatCapacityRatio / (heatCapacityRatio - 1.)); }
inline real calc_h(real rho, real P) { return calc_H(P) / rho; }
inline real calc_hTotal(real rho, real P, real ETotal) { return (P + ETotal) / rho; }
inline real calc_HTotal(real P, real ETotal) { return P + ETotal; }
inline real calc_eKin(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.v, x); }
inline real calc_EKin(<?=eqn.prim_t?> W, real3 x) { return W.rho * calc_eKin(W, x); }
inline real calc_EInt(<?=eqn.prim_t?> W) { return W.P / (heatCapacityRatio - 1.); }
inline real calc_eInt(<?=eqn.prim_t?> W) { return calc_EInt(W) / W.rho; }
inline real calc_EKin_fromCons(<?=eqn.cons_t?> U, real3 x) { return .5 * coordLenSq(U.m, x) / U.rho; }
inline real calc_ETotal(<?=eqn.prim_t?> W, real3 x) {
	real EPot = W.rho * W.ePot;
	return calc_EKin(W, x) + calc_EInt(W) + EPot;
}

inline real calc_Cs(const <?=eqn.prim_t?>* W) {
	return sqrt(heatCapacityRatio * W->P / W->rho);
}
]], {
		eqn = self,
	})
end

function Euler:getPrimConsCode()
	return template([[

inline <?=eqn.prim_t?> primFromCons(<?=eqn.cons_t?> U, real3 x) {
	real EPot = U.rho * U.ePot;
	real EKin = calc_EKin_fromCons(U, x);
	real EInt = U.ETotal - EKin - EPot;
	return (<?=eqn.prim_t?>){
		.rho = U.rho,
		.v = real3_real_mul(U.m, 1./U.rho),
		.P = (heatCapacityRatio - 1.) * EInt,
		.ePot = U.ePot,
	};
}

inline <?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> W, real3 x) {
	return (<?=eqn.cons_t?>){
		.rho = W.rho,
		.m = real3_real_mul(W.v, W.rho),
		.ETotal = calc_ETotal(W, x),
		.ePot = W.ePot,
	};
}

<?=eqn.cons_t?> apply_dU_dW(
	<?=eqn.prim_t?> WA, 
	<?=eqn.prim_t?> W, 
	real3 x
) {
	real3 WA_vL = coord_lower(WA.v, x);
	return (<?=eqn.cons_t?>){
		.rho = W.rho,
		.m = real3_add(
			real3_real_mul(WA.v, W.rho), 
			real3_real_mul(W.v, WA.rho)),
		.ETotal = W.rho * .5 * real3_dot(WA.v, WA_vL) 
			+ WA.rho * real3_dot(W.v, WA_vL)
			+ W.P / (heatCapacityRatio - 1.)
			+ WA.rho * W.ePot,
		.ePot = W.ePot,
	};
}

<?=eqn.prim_t?> apply_dW_dU(
	<?=eqn.prim_t?> WA,
	<?=eqn.cons_t?> U,
	real3 x
) {
	real3 WA_vL = coord_lower(WA.v, x);
	return (<?=eqn.prim_t?>){
		.rho = U.rho,
		.v = real3_sub(
			real3_real_mul(U.m, 1. / WA.rho),
			real3_real_mul(WA.v, U.rho / WA.rho)),
		.P = (heatCapacityRatio - 1.) * (
			.5 * real3_dot(WA.v, WA_vL) * U.rho 
			- real3_dot(U.m, WA_vL)
			+ U.ETotal 
			- WA.rho * U.ePot),
		.ePot = U.ePot,
	};
}


]], {
		eqn = self,
	})
end

Euler.initStateCode = [[
<? local xNames = require 'common'().xNames ?>
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_real_mul(real3_add(mins, maxs), .5);
	bool lhs = true
<?
for i=1,solver.dim do
	local xi = xNames[i]
?>	&& x.<?=xi?> < mids.<?=xi?>
<?
end
?>;
	
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	
	//TODO make this B for Maxwell
	
	real3 B = real3_zero;	//set for MHD / thrown away for pure Euler
	real ePot = 0;

	<?=code?>

	<?=eqn.prim_t?> W = {
		.rho = rho,
		.v = cartesianToCoord(v, x),	//transform from cartesian to coordinate space 
		.P = P,
		.ePot = ePot,
	};
	UBuf[index] = consFromPrim(W, x);
}
]]

Euler.solverCodeFile = 'eqn/euler.cl'

Euler.displayVarCodeUsesPrims = true

-- k is 0,1,2
local function vorticity(eqn,k,result)
	local xs = {'x','y','z'}
	local i = (k+1)%3
	local j = (i+1)%3
	return {['vorticity '..xs[k+1]] = template([[
	if (OOB(1,1)) {
		<?=result?> = 0.;
	} else {
		global const <?=eqn.cons_t?>* Uim = U - stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Uip = U + stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Ujm = U - stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + stepsize.s<?=j?>;

		//TODO incorporate metric

		real vim_j = Uim->m.s<?=j?> / Uim->rho;
		real vip_j = Uip->m.s<?=j?> / Uip->rho;
		
		real vjm_i = Ujm->m.s<?=i?> / Ujm->rho;
		real vjp_i = Ujp->m.s<?=i?> / Ujp->rho;
		
		<?=result?> = (vjp_i - vjm_i) / (2. * grid_dx<?=i?>)
				- (vip_j - vim_j) / (2. * grid_dx<?=j?>);
	}
]], {
		i = i,
		j = j,
		eqn = eqn,
		result = result,
	})}
end

function Euler:getDisplayVars()
	local vars = Euler.super.getDisplayVars(self)
	vars:append{
		{v = '*value_real3 = W.v;', type='real3'},
		{P = '*value = W.P;'},
		{eInt = '*value = calc_eInt(W);'},
		{eKin = '*value = calc_eKin(W, x);'},
		{eTotal = '*value = U->ETotal / W.rho;'},
		{EInt = '*value = calc_EInt(W);'},
		{EKin = '*value = calc_EKin(W, x);'},
		{EPot = '*value = U->rho * U->ePot;'},
		{S = '*value = W.P / pow(W.rho, (real)heatCapacityRatio);'},
		{H = '*value = calc_H(W.P);'},
		{h = '*value = calc_h(W.rho, W.P);'},
		{HTotal = '*value = calc_HTotal(W.P, U->ETotal);'},
		{hTotal = '*value = calc_hTotal(W.rho, W.P, U->ETotal);'},
		{['Speed of Sound'] = '*value = calc_Cs(&W);'},
		{['Mach number'] = '*value = coordLen(W.v, x) / calc_Cs(&W);'},
	}:append{self.gravOp and
		{gravity = template([[
	if (OOB(1,1)) {
		*value = 0.;
	} else {
		<? 
for side=0,solver.dim-1 do ?>{
			global const <?=eqn.cons_t?>* Um = U - stepsize.s<?=side?>;
			global const <?=eqn.cons_t?>* Up = U + stepsize.s<?=side?>;
			value_real3->s<?=side?> = -(Up-><?=eqn.gravOp.potentialField?> - Um-><?=eqn.gravOp.potentialField?>) / (2. * dx<?=side?>_at(i));
		}<? 
end
for side=solver.dim,2 do ?>
		value_real3->s<?=side?> = 0.;
<? end ?>
	}
]], {eqn=self, solver=self.solver}), type='real3'} or nil
	}:append{
		{['temp (C)'] = template([[
<? local clnumber = require 'cl.obj.number' ?>
<? local materials = require 'materials' ?>
#define _0_C_in_K		273.15
#define C_v				<?=('%.50f'):format(materials.Air.C_v)?>
	*value = calc_eInt(W) / C_v - _0_C_in_K;
]])}
	}

	-- vorticity = [,x ,y ,z] [v.x, v.y, v.z][
	-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
	
	if not require 'solver.meshsolver'.is(self.solver) then
		if self.solver.dim == 2 then
			vars:insert(vorticity(self,2,'*value'))
		elseif self.solver.dim == 3 then
			local v = range(0,2):map(function(i) return vorticity(self,i,'value['..i..']') end)
			vars:insert{vorticityVec = template([[
	<? for i=0,2 do ?>{
		<?=select(2,next(v[i+1]))?>
	}<? end ?>
]], {v=v}), type='real3'}
		end
	end

	return vars
end

Euler.eigenVars = table{
	-- Roe-averaged vars
	{rho = 'real'},
	{v = 'real3'},
	{hTotal = 'real'},
	-- derived vars
	{vSq = 'real'},
	{Cs = 'real'},
}

function Euler:eigenWaveCodePrefix(side, eig, x)
	return template([[
	real Cs_sqrt_gU = <?=eig?>.Cs * coord_sqrt_gU<?=side..side?>(<?=x?>);
	real v_n = <?=eig?>.v.s[<?=side?>];
]], {
		eig = '('..eig..')',
		side = side,
		x = x,
	})
end

-- W is an extra param specific to Euler's calcDT in this case
function Euler:consWaveCodePrefix(side, U, x, W)
	return template([[
<? if not W then ?>
	<?=eqn.prim_t?> W = primFromCons(<?=U?>, <?=x?>);
<? end ?>
	real Cs_sqrt_gU = calc_Cs(&<?=W or 'W'?>) * coord_sqrt_gU<?=side..side?>(<?=x?>);
	real v_n = <?=W or 'W'?>.v.s[<?=side?>];
]], {
		eqn = self,
		U = '('..U..')',
		W = W and '('..W..')' or nil,
		side = side,
		x = x,
	})
end

function Euler:consWaveCode(side, U, x, waveIndex)
	if waveIndex == 0 then
		return '(v_n - Cs_sqrt_gU)'
	elseif waveIndex >= 1 and waveIndex <= 3 then
		return 'v_n'
	elseif waveIndex == 4 then
		return '(v_n + Cs_sqrt_gU)'
	end
	error'got a bad waveIndex'
end

-- as long as U or eig isn't used, we can use this for both implementations
Euler.eigenWaveCode = Euler.consWaveCode

-- this one calcs cell prims once and uses it for all sides
-- it is put here instead of in eqn/euler.cl so euler-burgers can override it
function Euler:getCalcDTCode()
	return template([[
<? local solver = eqn.solver ?>
kernel void calcDT(
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 x = cell_x(i);

	const global <?=eqn.cons_t?>* U = UBuf + index;
	<?=eqn.prim_t?> W = primFromCons(*U, x);
	real Cs = calc_Cs(&W);

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		//use cell-centered eigenvalues
		real lambdaMin = W.v.s<?=side?> - Cs;
		real lambdaMax = W.v.s<?=side?> + Cs;
		lambdaMin = min((real)0., lambdaMin);
		lambdaMax = max((real)0., lambdaMax);
		dt = min(dt, ((real)grid_dx<?=side?> / (fabs(lambdaMax - lambdaMin) + (real)1e-9)));
	}<? end ?>
	dtBuf[index] = dt;
}
]], {
		eqn = self,
	})
end

return Euler

--[[
scratch work:

P V = n R T = m R_spec T
ideal gas EOS: P V_m = R T
T = T_C + 273.15' C
R = gas constant = 8.314459848 J / (mol K) 
R_spec = 287.058 J / (kg K) for dry air 
R_spec = R / M = k_B / m
M = molar mass
k_B = Boltzmann constant = 1.3806485279e-23 J / K
m = mass within volume V	<-> rho = m / V
R_spec = k_B rho V

caloric perfect:
P = (gamma - 1) rho e_int
gamma = C_p / C_v = adiabatic index / ratio of specific heats
e_int = C_v T = internal energy per unit mass / specific internal energy
C_v = specific heat at constant volume
C_p = specific heat at constant pressure
P = rho R_spec T
e_int = C_v T
P = rho (C_p - C_v) e_int / C_v = rho (C_p / C_v - 1) e_int

0 C = 273.15 K

using some real-world numbers ...
P = (gamma - 1) rho e_int
101325 kg / (m s^2) = (C_p - C_v) / C_v (1.2754 kg / m^3) C_v T 
101325 kg / (m s^2) = (1006 - 717.1) J / (kg K) (1.2754 kg / m^3) T 
T = 274.99364522457 K
T = 1.8436452245715 C ... should be


Source terms:

rho,t + (rho v^j)_,j = 0
(rho v^i),t + (rho v^i v^j + g^ij P)_,j = 0
E_total,t + (v^j H_total)_,j = 0

m^i = rho v^i
E_total = E_kin + E_int <=> E_int = E_total - E_kin
P = (gamma - 1) E_int = (gamma - 1) (E_total - E_kin) = (gamma - 1) (E_total - 1/2 m^2 / rho)
H_total = E_total + P = gamma E_total - 1/2 (gamma - 1) m^2 / rho

rho,t + m^j_,j = 0
m^i_,t + (m^i m^j / rho + g^ij (gamma - 1) (E_total - 1/2 m^2 / rho))_,j = 0
E_total,t + (m^j (gamma E_total / rho - 1/2 (gamma - 1) m^2 / rho^2))_,j = 0

rho,t + m^j_,j = 0
m^i_,t + m^i_,j m^j / rho + m^i m^j_,j / rho - m^i m^j rho_,j / rho^2 
	+ g^ij_,j (gamma - 1) (E_total - 1/2 m^2 / rho) 
	+ g^ij (gamma - 1) (E_total_,j - 1/2 m^2 / rho) 
	+ g^ij (gamma - 1) (E_total - m^k_,j m_k / rho) 
	+ g^ij (gamma - 1) (E_total - 1/2 m^k m^l g_kl,j / rho) 
	+ g^ij (gamma - 1) (E_total + 1/2 m^2 / rho^2 rho_,j)
	= 0
E_total,t + m^j_,j (gamma E_total / rho - 1/2 (gamma - 1) m^2 / rho^2) 
	+ m^j (gamma E_total_,j / rho - 1/2 (gamma - 1) m^2 / rho^2) 
	+ m^j (-gamma E_total / rho^2 rho_,j - 1/2 (gamma - 1) m^2 / rho^2) 
	+ m^j (gamma E_total / rho - (gamma - 1) m^k_,j m_k / rho^2) 
	+ m^j (gamma E_total / rho - 1/2 (gamma - 1) m^k m^l g_kl,j / rho^2) 
	+ m^j (gamma E_total / rho + (gamma - 1) m^2 rho_,j / rho^3) 
	= 0

rho,t + m^j_,j = 0
m^i_,t 
	+ ((1/2 (gamma - 1) g^ij m^2 - m^i m^j) / rho^2) rho_,j
	+ ((delta^i_k m^j + m^i delta^j_k - (gamma - 1) g^ij m_k) / rho) m^k_,j 
	+ g^ij (gamma - 1) E_total_,j
	+ (- g^ik g^lj P - 1/2 (gamma - 1) g^ij m^k m^l / rho) g_kl,j
	 = 0
E_total,t 
	+ (( (gamma - 1) m^j m^2 / rho - gamma E_total) / rho^2) rho_,j
	+ (delta^j_k h_total - (gamma - 1) m^j m_k / rho^2) m^k_,j
	+ (gamma m^j / rho) E_total_,j 
	+ (-1/2 (gamma - 1) m^j m^k m^l / rho^2) g_kl,j
	 = 0

rho,t + m^j_,j = 0
m^i_,t 
	+ (1/2 (gamma - 1) g^ij v^2 - v^i v^j) rho_,j
	+ (delta^i_k v^j + delta^j_k v^i - (gamma - 1) g^ij v_k) m^k_,j 
	+ g^ij (gamma - 1) E_total_,j
	+ (-g^ik g^lj P - 1/2 (gamma - 1) g^ij rho v^k v^l) g_kl,j
	 = 0
E_total,t 
	+ ( (gamma - 1) v^j v^2 - gamma E_total / rho^2) rho_,j
	+ (delta^j_k h_total - (gamma - 1) v^j v_k) m^k_,j
	+ gamma v^j E_total_,j
	+ (-1/2 (gamma - 1) rho v^j v^k v^l) g_kl,j
	 = 0

rho,t + m^j_,j + Gamma^j_kj m^k = 0
m^i_,t 
	+ (1/2 (gamma - 1) g^ij v^2 - v^i v^j) rho_,j
	+ (delta^i_k v^j + delta^j_k v^i - (gamma - 1) g^ij v_k) (m^k_,j + Gamma^k_lkj m^l)
	+ g^ij (gamma - 1) E_total_,j
	 = 0
E_total,t 
	+ ( (gamma - 1) v^j v^2 - gamma E_total / rho^2) rho_,j
	+ (delta^j_k h_total - (gamma - 1) v^j v_k) (m^k_,j + Gamma^k_lj m^l)
	+ gamma v^j E_total_,j
	 = 0

U^I_,t + F^Ij_,j + ...
	... + Gamma^j_kj m^k = 0
	... + (delta^i_k v^j + delta^j_k v^i - (gamma - 1) g^ij v_k) Gamma^k_lj m^l + (g^ik g^lj P + 1/2 (gamma - 1) g^ij m^k m^l / rho) (Gamma_klj + Gamma_lkj) = 0
	... + (delta^j_k h_total - (gamma - 1) v^j v_k) Gamma^k_lj m^l + (1/2 (gamma - 1) m^j m^k m^l / rho^2) (Gamma_klj + Gamma_lkj) = 0

U^I_,t + 1/sqrt(g) ( sqrt(g) F^Ij )_,j + ...
	...   [                  0                 ]
	... = [ -Gamma^i_jk (rho v^j v^k + P g^jk) ]
	...   [                  0                 ]


~

Or we could remove the metric from the flux:

[  rho  ]     [1   0  0] [          rho v^j        ]      [ 0 ]
[rho v^i]   + [0 g^ik 0] [rho v_k v^j + delta_k^j P]    = [ 0 ]
[E_total],t   [0   0  1] [        v^j H_total      ]_,j   [ 0 ]

the metric can be factored outside of the time partial derivatives, but only for static grids

now we have a delta_k^j instead of a g^ij, but at the cost that the v_k v^j is no longer symmetric as the v^i v^j term would be
--]]

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

Euler.hasEigenCode = true
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
	self.gravOp = SelfGrav{solver = self.solver}
	self.solver.ops:insert(self.gravOp)
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
--in order to make things work, gamma needs to be set *HERE AND IN INIT/EULER*
-- which means it is being read and written in multiple places
-- TODO consolidate that
	self:addGuiVars{
		{name='heatCapacityRatio', value=7/5},
		{name='useNavierStokesViscosityTerm', value=false},
		{name='viscosity_K', value=.1},
		{name='viscosity_mu_T', value=1e-3},
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
		.v = real3_scale(U.m, 1./U.rho),
		.P = (heatCapacityRatio - 1.) * EInt,
		.ePot = U.ePot,
	};
}

inline <?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> W, real3 x) {
	return (<?=eqn.cons_t?>){
		.rho = W.rho,
		.m = real3_scale(W.v, W.rho),
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
			real3_scale(WA.v, W.rho), 
			real3_scale(W.v, WA.rho)),
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
			real3_scale(U.m, 1. / WA.rho),
			real3_scale(WA.v, U.rho / WA.rho)),
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
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	bool lhs = true
<?
for i=1,solver.dim do
	local xi = xNames[i]
?>	&& x.<?=xi?> < mids.<?=xi?>
<?
end
?>;
	
	real rho = 0;
	real3 v = _real3(0,0,0);
	real P = 0;
	
	//TODO make this B for Maxwell
	
	real3 B = _real3(0,0,0);	//set for MHD / thrown away for pure Euler
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
		{v = '*valuevec = W.v;', type='real3'},
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
			valuevec->s<?=side?> = -(Up-><?=eqn.gravOp.potentialField?> - Um-><?=eqn.gravOp.potentialField?>) / (2. * dx<?=side?>_at(i));
		}<? 
end
for side=solver.dim,2 do ?>
		valuevec->s<?=side?> = 0.;
<? end ?>
	}
]], {eqn=self, solver=self.solver}), type='real3'} or nil}

	-- vorticity = [,x ,y ,z] [v.x, v.y, v.z][
	-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
		
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
	real Cs_sqrt_gU = <?=eig?>->Cs * coord_sqrt_gU<?=side..side?>(<?=x?>);
	real v_n = <?=eig?>->v.s[<?=side?>];
]], {
		eig = '('..eig..')',
		side = side,
		x = x,
	})
end

function Euler:eigenWaveCode(side, eig, x, waveIndex)
	if waveIndex == 0 then
		return '(v_n - Cs_sqrt_gU)'
	elseif waveIndex >= 1 and waveIndex <= 3 then
		return 'v_n'
	elseif waveIndex == 4 then
		return '(v_n + Cs_sqrt_gU)'
	end
	error'got a bad waveIndex'
end

function Euler:getFluxFromConsCode()
	return template([[
<? local solver = eqn.solver ?>
<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(U, x);
	real vj = W.v.s<?=side?>;
	real HTotal = U.ETotal + W.P;
	
	<?=eqn.cons_t?> F;
	F.rho = U.m.s<?=side?>;
	F.m = real3_scale(U.m, vj);
<? for i=0,2 do
?>	F.m.s<?=i?> += coord_gU<?=i?><?=side?>(x) * W.P;
<? end
?>	F.ETotal = HTotal * vj;
	F.ePot = 0;
	return F;
}
<? end ?>
]], {
		eqn = self,
	})
end

return Euler

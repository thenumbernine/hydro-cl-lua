--[[
Taken from Petrova - Finite Volume Methods: Powerful Means of Engineering Design

k-omega turbulence model of Navier Stokes method for finite volume
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local template = require 'template'
local materials = require 'materials'
local Equation = require 'eqn.eqn'

local NavierStokesWilcox = class(Equation)
NavierStokesWilcox.name = 'Navier-Stokes -- 1998 Wilcox'

NavierStokesWilcox.numStates = 8	-- ePot is the last param -- specific potential energy
NavierStokesWilcox.numWaves = 7	-- v-a, v,v,v,v,v, v+a
NavierStokesWilcox.numIntStates = 7

NavierStokesWilcox.hasEigenCode = true
NavierStokesWilcox.hasFluxFromConsCode = true
NavierStokesWilcox.roeUseFluxFromCons = true
NavierStokesWilcox.useSourceTerm = true

NavierStokesWilcox.initStates = require 'init.euler'

function NavierStokesWilcox:init(args)
	NavierStokesWilcox.super.init(self, args)

--[[ might have to change the selfgrav terms ...
	local SelfGrav = require 'op.selfgrav'

if require 'solver.meshsolver'.is(self.solver) then
	print("not using selfgrav with mesh solvers yet")
else
		self.gravOp = SelfGrav{solver = self.solver}
		self.solver.ops:insert(self.gravOp)
end
--]]
end

NavierStokesWilcox.primVars = table{
	{name='rhoBar', type='real'},
	{name='vTilde', type='real3'},
	{name='PStar', type='real'},
	{name='k', type='real'},
	{name='omega', type='real'},
	{name='ePot', type='real'},
}

NavierStokesWilcox.consVars = table{
	{name='rhoBar', type='real'},
	{name='rhoBar_vTilde', type='real3'},
	{name='rhoBar_eTotalTilde', type='real'},
	{name='rhoBar_k', type='real'},
	{name='rhoBar_omega', type='real'},
	{name='ePot', type='real'},
}

function NavierStokesWilcox:createInitState()
	NavierStokesWilcox.super.createInitState(self)
--in order to make things work, gamma needs to be set *HERE AND IN INIT/EULER*
-- which means it is being read and written in multiple places
-- TODO consolidate that
	self:addGuiVars{
		{name='C_v', value=materials.Air.C_v},
		{name='C_p', value=materials.Air.C_p},

		-- specific gas constant R_spec
		{name='gasConstant', value=(materials.Air.C_p - materials.Air.C_v)},
		
		{name='heatCapacityRatio', value=materials.Air.C_p / materials.Air.C_v},
	}
end

function NavierStokesWilcox:getCommonFuncCode()
	return template([[
#define R_over_C_v (solver->gasConstant / solver->C_v)
#define C_v_over_R (solver->C_v / solver->gasConstant)

//inline real calc_H(real PStar) { return PStar * ((R_over_C_v + 1.) / (R_over_C_v)); }
//inline real calc_h(real rhoBar, real PStar) { return calc_H(PStar) / rhoBar; }
//inline real calc_hTotal(real rhoBar, real PStar, real rhoBar_eTotalTilde) { return (PStar + rhoBar_eTotalTilde) / rhoBar; }
//inline real calc_HTotal(real PStar, real rhoBar_eTotalTilde) { return PStar + rhoBar_eTotalTilde; }
static inline real calc_eKinTilde(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.vTilde, x); }
static inline real calc_EKinTilde(<?=eqn.prim_t?> W, real3 x) { return W.rhoBar * calc_eKinTilde(W, x); }

//before
//inline real calc_EIntTilde(<?=eqn.prim_t?> W) { return W.PStar * C_v_over_R; }
//inline real calc_eIntTilde(<?=eqn.prim_t?> W) { return calc_EIntTilde(W) / W.rhoBar; }

//after
static inline real calc_PBar(<?=eqn.prim_t?> W) { return W.PStar - 2./3. * W.rhoBar * W.k; }
static inline real calc_TTilde(<?=eqn.prim_t?> W) { return calc_PBar(W) / (W.rhoBar * solver->gasConstant); }
static inline real calc_eIntTilde(<?=eqn.prim_t?> W) { return solver->C_v * calc_TTilde(W); }
static inline real calc_EIntTilde(<?=eqn.prim_t?> W) { return W.rhoBar * calc_eIntTilde(W); }

static inline real calc_EKin_fromCons(<?=eqn.cons_t?> U, real3 x) { return .5 * coordLenSq(U.rhoBar_vTilde, x) / U.rhoBar; }
static inline real calc_ETotal(<?=eqn.prim_t?> W, real3 x) {
	return calc_EKinTilde(W, x) + calc_EIntTilde(W);
}

static inline real calc_Cs(const <?=eqn.prim_t?> W) {
	return sqrt((R_over_C_v + 1.) * W.PStar / W.rhoBar);
}
]], {
		eqn = self,
	})
end

function NavierStokesWilcox:getPrimConsCode()
	return template([[
static inline <?=eqn.prim_t?> primFromCons(<?=eqn.cons_t?> U, real3 x) {
	real3 vTilde = real3_real_mul(U.rhoBar_vTilde, 1. / U.rhoBar);
	real vTildeSq = coordLenSq(vTilde, x);
	real rhoBar_eIntTilde = U.rhoBar_eTotalTilde - .5 * U.rhoBar * vTildeSq - U.rhoBar_k;
	real rhoBar_TTilde = rhoBar_eIntTilde / solver->C_v;
	real PBar = rhoBar_TTilde * solver->gasConstant;
	real PStar = PBar + 2./3. * U.rhoBar_k;
	
	return (<?=eqn.prim_t?>){
		.rhoBar = U.rhoBar,
		.vTilde = vTilde,
		.PStar = PStar,
		.k = U.rhoBar_k / U.rhoBar,
		.omega = U.rhoBar_omega / U.rhoBar,
		.ePot = U.ePot,
	};
}

static inline <?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> W, real3 x) {
	real rhoBar_k = W.rhoBar * W.k;

	//eqn 6: PStar = PBar + 2/3 rhoBar k
	real PBar = W.PStar - 2./3. * rhoBar_k;

	//eqn 10: PBar = rhoBar R TTilde
	real TTilde = PBar / (W.rhoBar * solver->gasConstant);
	
	//eqn 6: eIntTilde = C_v TTilde
	real eIntTilde = solver->C_v * TTilde;
	
	//eqn 6: eTotalTilde = eIntTilde + 1/2 vTilde^2 + W.k
	//so eTotalTilde = C_v PStar / rhoBar + 1/2 vTilde^2 + (1 - 2/3 C_v / solver->gasConstant) k
	real eTotalTilde = eIntTilde + .5 * coordLenSq(W.vTilde, x) + W.k;
	
	return (<?=eqn.cons_t?>){
		.rhoBar = W.rhoBar,
		.rhoBar_vTilde = real3_real_mul(W.vTilde, W.rhoBar),
		.rhoBar_eTotalTilde = W.rhoBar * eTotalTilde,
		.rhoBar_k = rhoBar_k,
		.rhoBar_omega = W.rhoBar * W.omega,
		.ePot = W.ePot,
	};
}

<?=eqn.cons_t?> apply_dU_dW(
	<?=eqn.prim_t?> WA, 
	<?=eqn.prim_t?> W, 
	real3 x
) {
	real3 WA_vTildeL = coord_lower(WA.vTilde, x);
	return (<?=eqn.cons_t?>){
		.rhoBar = W.rhoBar,
		.rhoBar_vTilde = real3_add(
			real3_real_mul(WA.vTilde, W.rhoBar), 
			real3_real_mul(W.vTilde, WA.rhoBar)),
		.rhoBar_eTotalTilde = W.rhoBar * (.5 * real3_dot(WA.vTilde, WA_vTildeL) 
				+ (1. - 2./3. * C_v_over_R) * WA.k)
			+ WA.rhoBar * real3_dot(W.vTilde, WA_vTildeL)
			+ W.PStar * C_v_over_R
			+ (1. - 2./3 * C_v_over_R) * WA.rhoBar * W.k,
		.rhoBar_k = WA.k * W.rhoBar + WA.rhoBar * W.k,
		.rhoBar_omega = WA.omega * W.rhoBar + WA.rhoBar * W.omega,
		.ePot = W.ePot,
	};
}

<?=eqn.prim_t?> apply_dW_dU(
	<?=eqn.prim_t?> WA,
	<?=eqn.cons_t?> U,
	real3 x
) {
	real3 WA_vTildeL = coord_lower(WA.vTilde, x);
	return (<?=eqn.prim_t?>){
		.rhoBar = U.rhoBar,
		.vTilde = real3_sub(
			real3_real_mul(U.rhoBar_vTilde, 1. / WA.rhoBar),
			real3_real_mul(WA.vTilde, U.rhoBar / WA.rhoBar)),
		.PStar = R_over_C_v * (
				.5 * real3_dot(WA.vTilde, WA_vTildeL) * U.rhoBar 
				- real3_dot(U.rhoBar_vTilde, WA_vTildeL)
				+ U.rhoBar_eTotalTilde
			) + (2./3. * R_over_C_v - 1.) * U.rhoBar_k,
		.k = U.rhoBar_k / WA.rhoBar - WA.k / WA.rhoBar * U.rhoBar,
		.omega = U.rhoBar_omega / WA.rhoBar - WA.omega / WA.rhoBar * U.rhoBar,
		.ePot = U.ePot,
	};
}


]], {
		eqn = self,
	})
end

NavierStokesWilcox.initStateCode = [[
<?
local common = require 'common'
local xNames = common.xNames
?>
kernel void initState(
	constant <?=solver.solver_t?>* solver,
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
	
	real3 B = real3_zero;	//set for MHD / thrown away for pure NavierStokesWilcox
	real ePot = 0;

	<?=code?>

	<?=eqn.prim_t?> W = {
		.rhoBar = rho,
		.vTilde = cartesianToCoord(v, x),	//transform from cartesian to coordinate space 
		.PStar = P,
		.k = 0,
		.omega = 0,
		.ePot = ePot,
	};
	UBuf[index] = consFromPrim(W, x);
}
]]

NavierStokesWilcox.solverCodeFile = 'eqn/navstokes-wilcox.cl'

NavierStokesWilcox.displayVarCodeUsesPrims = true

-- k is 0,1,2
local function vorticity(eqn,k,result)
	local xs = {'x','y','z'}
	local i = (k+1)%3
	local j = (i+1)%3
	return {
		name = 'vorticity '..xs[k+1], 
		code = template([[
	if (OOB(1,1)) {
		<?=result?> = 0.;
	} else {
		global const <?=eqn.cons_t?>* Uim = U - solver->stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Uip = U + solver->stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Ujm = U - solver->stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + solver->stepsize.s<?=j?>;
		
		//TODO incorporate metric
		
		real vim_j = Uim->rhoBar_vTilde.s<?=j?> / Uim->rhoBar;
		real vip_j = Uip->rhoBar_vTilde.s<?=j?> / Uip->rhoBar;
		
		real vjm_i = Ujm->rhoBar_vTilde.s<?=i?> / Ujm->rhoBar;
		real vjp_i = Ujp->rhoBar_vTilde.s<?=i?> / Ujp->rhoBar;
		
		<?=result?> = (vjp_i - vjm_i) / (2. * solver->grid_dx.s<?=i?>)
				- (vip_j - vim_j) / (2. * solver->grid_dx.s<?=j?>);
	}
]], {
		i = i,
		j = j,
		eqn = eqn,
		result = result,
	})}
end

function NavierStokesWilcox:getDisplayVars()
	local vars = NavierStokesWilcox.super.getDisplayVars(self)
	vars:append{
		{name='vTilde', code='value.vreal3 = W.vTilde;', type='real3'},
		{name='PStar', code='value.vreal = W.PStar;'},
		{name='eIntTilde', code='value.vreal = calc_eIntTilde(W);'},
		{name='eKinTilde', code='value.vreal = calc_eKinTilde(W, x);'},
		{name='eTotalTilde', code='value.vreal = U->rhoBar_eTotalTilde / W.rhoBar;'},
		{name='EIntTilde', code='value.vreal = calc_EIntTilde(W);'},
		{name='EKinTilde', code='value.vreal = calc_EKinTilde(W, x);'},
		{name='EPot', code='value.vreal = U->rhoBar * U->ePot;'},
		{name='S', code='value.vreal = W.PStar / pow(W.rhoBar, R_over_C_v + 1. );'},
		--{name='H', code='value.vreal = calc_H(W.PStar);'},
		--{name='h', code='value.vreal = calc_h(W.rhoBar, W.PStar);'},
		--{name='HTotal', code='value.vreal = calc_HTotal(W.PStar, U->rhoBar_eTotalTilde);'},
		--{name='hTotal', code='value.vreal = calc_hTotal(W.rhoBar, W.PStar, U->rhoBar_eTotalTilde);'},
		{name='Speed of Sound', code='value.vreal = calc_Cs(W);'},
		--{name='Mach number', code='value.vreal = coordLen(W.vTilde, x) / calc_Cs(W);'},
	}:append{self.gravOp and
		{name='gravity', code=template([[
	if (OOB(1,1)) {
		value.vreal = 0.;
	} else {
		<? 
for side=0,solver.dim-1 do ?>{
			global const <?=eqn.cons_t?>* Um = U - solver->stepsize.s<?=side?>;
			global const <?=eqn.cons_t?>* Up = U + solver->stepsize.s<?=side?>;
			value_real3->s<?=side?> = -(Up-><?=eqn.gravOp.potentialField?> - Um-><?=eqn.gravOp.potentialField?>) / (2. * cell_dx<?=side?>(x));
		}<? 
end
for side=solver.dim,2 do ?>
		value_real3->s<?=side?> = 0.;
<? end ?>
	}
]], {eqn=self, solver=self.solver}), type='real3'} or nil
	}:append{
		{name='temp', code='value.vreal = calc_eIntTilde(W) / solver->C_v;'},
	}

	-- vorticity = [,x ,y ,z] [vTilde.x, vTilde.y, vTilde.z][
	-- = [vTilde.z,y - vTilde.y,z; vTilde.x,z - vTilde.z,x; vTilde.y,x - vTilde.x,y]
	
	if not require 'solver.meshsolver'.is(self.solver) then
		if self.solver.dim == 2 then
			vars:insert(vorticity(self,2,'value.vreal'))
		elseif self.solver.dim == 3 then
			local vTilde = range(0,2):map(function(i) return vorticity(self,i,'value['..i..']') end)
			vars:insert{name='vorticityVec', code=template([[
	<? for i=0,2 do ?>{
		<?=select(2,next(vTilde[i+1]))?>
	}<? end ?>
]], {vTilde=vTilde}), type='real3'}
		end
	end

	return vars
end

NavierStokesWilcox.eigenVars = table{
	-- Roe-averaged vars
	{name='rhoBar', type='real'},
	{name='vTilde', type='real3'},
	{name='hTotal', type='real'},
	{name='k', type='real'},
	{name='omega', type='real'},
	-- derived vars
	{name='vTildeSq', type='real'},
	{name='Cs', type='real'},
}

function NavierStokesWilcox:eigenWaveCodePrefix(side, eig, x)
	return template([[
	real Cs_sqrt_gU = <?=eig?>->Cs * coord_sqrt_g_uu<?=side..side?>(<?=x?>);
	real v_n = <?=eig?>->vTilde.s[<?=side?>];
]], {
		eig = '('..eig..')',
		side = side,
		x = x,
	})
end

function NavierStokesWilcox:eigenWaveCode(side, eig, x, waveIndex)
	if waveIndex == 0 then
		return '(v_n - Cs_sqrt_gU)'
	elseif waveIndex >= 1 and waveIndex <= 5 then
		return 'v_n'
	elseif waveIndex == 6 then
		return '(v_n + Cs_sqrt_gU)'
	end
	error'got a bad waveIndex'
end

return NavierStokesWilcox

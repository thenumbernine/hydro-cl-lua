--[[
Taken from Petrova - Finite Volume Methods: Powerful Means of Engineering Design

k-omega turbulence model of Navier Stokes method for finite volume
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local materials = require 'hydro.materials'
local Equation = require 'hydro.eqn.eqn'

local NavierStokesWilcox = class(Equation)
NavierStokesWilcox.name = 'Navier-Stokes -- 1998 Wilcox'

NavierStokesWilcox.numWaves = 7	-- v-a, v,v,v,v,v, v+a
NavierStokesWilcox.numIntStates = 7

NavierStokesWilcox.roeUseFluxFromCons = true
NavierStokesWilcox.useSourceTerm = true

NavierStokesWilcox.initConds = require 'hydro.init.euler'

function NavierStokesWilcox:init(args)
	self.primVars = table{
		{name='rhoBar', type='real'},
		{name='vTilde', type='real3'},
		{name='PStar', type='real'},
		{name='k', type='real'},
		{name='omega', type='real'},
		{name='ePot', type='real'},
	}

	self.consVars = table{
		{name='rhoBar', type='real'},
		{name='rhoBar_vTilde', type='real3'},
		{name='rhoBar_eTotalTilde', type='real'},
		{name='rhoBar_k', type='real'},
		{name='rhoBar_omega', type='real'},
		{name='ePot', type='real'},
	}

	if args.incompressible then
		self.consVars:insert{name='mPot', type='real', units='kg/(m*s)'}
		self.primVars:insert{name='mPot', type='real', units='kg/(m*s)'}
	end

	NavierStokesWilcox.super.init(self, args)


	if require 'hydro.solver.meshsolver'.is(self.solver) then
		print("not using selfgrav with mesh solvers yet")
	else
--[[ might have to change the selfgrav terms ...
		local SelfGrav = require 'hydro.op.selfgrav'
		self.gravOp = SelfGrav{solver = self.solver}
		self.solver.ops:insert(self.gravOp)
--]]
		if args.incompressible then
			local NoDiv = require 'hydro.op.nodiv'{
				poissonSolver = require 'hydro.op.poisson_jacobi',	-- krylov is having errors.  TODO bug in its boundary code?
			}
			self.solver.ops:insert(NoDiv{
				solver = self.solver,
				vectorField = 'rhoBar_vTilde',
				potentialField = 'mPot',
			
				-- div v = 0
				-- div (m/ρ) = 0
				-- 1/ρ div m - 1/ρ^2 m dot grad ρ = 0
				-- div m = (m dot grad ρ)/ρ 
				chargeCode = self:template[[
	<? for j=0,solver.dim-1 do ?>{
		global const <?=eqn.cons_t?>* Ujm = U - solver->stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + solver->stepsize.s<?=j?>;
		real drho_dx = (Ujp->rhoBar - Ujm->rhoBar) * (.5 / solver->grid_dx.s<?=j?>);
		source -= drho_dx * U->rhoBar_vTilde.s<?=j?> / U->rhoBar;
	}<? end ?>
]],
			})
		end
	end
end

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

function NavierStokesWilcox:initCodeModuleCommon()
	self.solver.modules:add{
		name = 'eqn.common',
		code = self:template[[
#define R_over_C_v (solver->gasConstant / solver->C_v)
#define C_v_over_R (solver->C_v / solver->gasConstant)

//real calc_H(real PStar) { return PStar * ((R_over_C_v + 1.) / (R_over_C_v)); }
//real calc_h(real rhoBar, real PStar) { return calc_H(PStar) / rhoBar; }
//real calc_hTotal(real rhoBar, real PStar, real rhoBar_eTotalTilde) { return (PStar + rhoBar_eTotalTilde) / rhoBar; }
//real calc_HTotal(real PStar, real rhoBar_eTotalTilde) { return PStar + rhoBar_eTotalTilde; }
real calc_eKinTilde(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.vTilde, x); }
real calc_EKinTilde(<?=eqn.prim_t?> W, real3 x) { return W.rhoBar * calc_eKinTilde(W, x); }

//before
//real calc_EIntTilde(<?=eqn.prim_t?> W) { return W.PStar * C_v_over_R; }
//real calc_eIntTilde(<?=eqn.prim_t?> W) { return calc_EIntTilde(W) / W.rhoBar; }

//after
real calc_PBar(<?=eqn.prim_t?> W) { return W.PStar - 2./3. * W.rhoBar * W.k; }
real calc_TTilde(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_PBar(W) / (W.rhoBar * solver->gasConstant); }
real calc_eIntTilde(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return solver->C_v * calc_TTilde(solver, W); }
real calc_EIntTilde(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return W.rhoBar * calc_eIntTilde(solver, W); }

real calc_EKin_fromCons(<?=eqn.cons_t?> U, real3 x) { return .5 * coordLenSq(U.rhoBar_vTilde, x) / U.rhoBar; }
real calc_ETotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) {
	return calc_EKinTilde(W, x) + calc_EIntTilde(solver, W);
}

real calc_Cs(constant <?=solver.solver_t?>* solver, const <?=eqn.prim_t?> W) {
	return sqrt((R_over_C_v + 1.) * W.PStar / W.rhoBar);
}
]],
	}
end

function NavierStokesWilcox:initCodeModulePrimCons()
	self.solver.modules:add{
		name = 'eqn.prim-cons',
		depends = {
			'eqn.prim_t',
			'eqn.cons_t',
			'coordLenSq',
		},
		code = self:template[[
<?=eqn.prim_t?> primFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x
) {
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

<?=eqn.cons_t?> consFromPrim(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> W,
	real3 x
) {
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
]],
	}

	-- only used by PLM
	self.solver.modules:add{
		name = 'eqn.dU-dW',
		depends = {
			'real3',
			'solver.solver_t',
			'eqn.prim_t',
			'eqn.cons_t',
		},
		code = self:template[[
<?=eqn.cons_t?> apply_dU_dW(
	constant <?=solver.solver_t?>* solver,
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
	constant <?=solver.solver_t?>* solver,
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
]],
	}
end

NavierStokesWilcox.initCondCode = [[
kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
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
	UBuf[index] = consFromPrim(solver, W, x);
}
]]

NavierStokesWilcox.solverCodeFile = 'hydro/eqn/navstokes-wilcox.cl'

function NavierStokesWilcox:getModuleDependsApplyInitCond() 
	return table(NavierStokesWilcox.super.getModuleDependsApplyInitCond(self)):append{
		'cartesianToCoord',
		'eqn.prim-cons',
	}
end

function NavierStokesWilcox:getModuleDependsSolver() 
	return table(NavierStokesWilcox.super.getModuleDependsSolver(self)):append{
		'eqn.common',
		'eqn.prim-cons',
		'coord_g_uu##',
		'coord_g_uu',
		'coord_sqrt_g_uu##',
		'coord_lower',
	}
end

NavierStokesWilcox.displayVarCodeUsesPrims = true

function NavierStokesWilcox:getDisplayVars()
	local vars = NavierStokesWilcox.super.getDisplayVars(self)
	vars:append{
		{name='vTilde', code='value.vreal3 = W.vTilde;', type='real3'},
		{name='PStar', code='value.vreal = W.PStar;'},
		{name='eIntTilde', code='value.vreal = calc_eIntTilde(solver, W);'},
		{name='eKinTilde', code='value.vreal = calc_eKinTilde(W, x);'},
		{name='eTotalTilde', code='value.vreal = U->rhoBar_eTotalTilde / W.rhoBar;'},
		{name='EIntTilde', code='value.vreal = calc_EIntTilde(solver, W);'},
		{name='EKinTilde', code='value.vreal = calc_EKinTilde(W, x);'},
		{name='EPot', code='value.vreal = U->rhoBar * U->ePot;'},
		{name='S', code='value.vreal = W.PStar / pow(W.rhoBar, R_over_C_v + 1. );'},
		--{name='H', code='value.vreal = calc_H(W.PStar);'},
		--{name='h', code='value.vreal = calc_h(W.rhoBar, W.PStar);'},
		--{name='HTotal', code='value.vreal = calc_HTotal(W.PStar, U->rhoBar_eTotalTilde);'},
		--{name='hTotal', code='value.vreal = calc_hTotal(W.rhoBar, W.PStar, U->rhoBar_eTotalTilde);'},
		{name='Speed of Sound', code='value.vreal = calc_Cs(solver, W);'},
		--{name='Mach number', code='value.vreal = coordLen(W.vTilde, x) / calc_Cs(solver, W);'},
	}:append{self.gravOp and
		{name='gravity', code=self:template[[
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
]], type='real3'} or nil
	}:append{
		{name='temp', code='value.vreal = calc_eIntTilde(solver, W) / solver->C_v;'},
	}

	vars:insert(self:createDivDisplayVar{
		field = 'vTilde', 
		getField = function(U, j)
			return U..'->rhoBar_vTilde.s'..j..' / '..U..'->rhoBar'
		end,
		units = 'kg/(m^3*s)',
	})

	vars:insert(self:createCurlDisplayVar{
		field = 'vTilde',
		getField = function(U, j)
			return U..'->rhoBar_vTilde.s'..j..' / '..U..'->rhoBar'
		end,
		units = 'm/s^2',
	})


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

function NavierStokesWilcox:eigenWaveCodePrefix(n, eig, x)
	return self:template([[
	real Cs_nLen = <?=eig?>.Cs * normal_len(<?=n?>);
	real v_n = normal_vecDotN1(<?=n?>, <?=eig?>.vTilde);
]], {
		eig = '('..eig..')',
		n = n,
		x = x,
	})
end

function NavierStokesWilcox:consWaveCodePrefix(n, U, x)
	return self:template([[
	<?=eqn.prim_t?> W = primFromCons(solver, <?=U?>, <?=x?>);
	real Cs_nLen = calc_Cs(solver, W) * normal_len(<?=n?>);
	real v_n = normal_vecDotN1(<?=n?>, W.vTilde);
]], {
		U = '('..U..')',
		n = n,
		x = x,
	})
end

function NavierStokesWilcox:consWaveCode(n, eig, x, waveIndex)
	if waveIndex == 0 then
		return '(v_n - Cs_nLen)'
	elseif waveIndex >= 1 and waveIndex <= 5 then
		return 'v_n'
	elseif waveIndex == 6 then
		return '(v_n + Cs_nLen)'
	end
	error'got a bad waveIndex'
end

NavierStokesWilcox.eigenWaveCode = NavierStokesWilcox.consWaveCode

return NavierStokesWilcox

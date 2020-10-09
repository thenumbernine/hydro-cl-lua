--[[
based on whatever my numerical_relativity_codegen z4 is based on, which is probably a Bona-Masso paper,
probably 2004 Bona et al "A symmetry-breaking mechanism for the Z4 general-covariant evolution system"
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local symmath = require 'symmath'
local EinsteinEqn = require 'hydro.eqn.einstein'
local Struct = require 'hydro.code.struct'

local common = require 'hydro.common'
local xNames = common.xNames
local symNames = common.symNames
local sym = common.sym


local Z4_2004Bona = class(EinsteinEqn)
Z4_2004Bona.name = 'Z4_2004Bona'
Z4_2004Bona.useSourceTerm = true
Z4_2004Bona.useConstrainU = true

--[[
args:
noZeroRowsInFlux = true by default.
	true = use 13 vars of a_x, d_xij, K_ij
	false = use all 30 or so hyperbolic conservation variables
useShift
	
	useShift = 'none'
	
	useShift = 'MinimalDistortionElliptic' -- minimal distortion elliptic via Poisson relaxation.  Alcubierre's book, eqn 4.3.14 and 4.3.15
	
	useShift = 'MinimalDistortionEllipticEvolve' -- minimal distortion elliptic via evolution.  eqn 10 of 1996 Balakrishna et al "Coordinate Conditions and their Implementations in 3D Numerical Relativity" 

	useShift = '2005 Bona / 2008 Yano'
	-- 2008 Yano et al, from 2005 Bona et al "Geometrically Motivated..."
	-- 2005 Bona mentions a few, but 2008 Yano picks the first one from the 2005 Bona paper.

	useShift = 'HarmonicShiftCondition-FiniteDifference'
	-- 2008 Alcubierre 4.3.37
	-- I see some problems in the warp bubble test ...

	useShift = 'LagrangianCoordinates'
	--[=[
	Step backwards along shift vector and advect the state
	Idk how accurate this is ...
	Hmm, even if I implement the Lie derivative as Lagrangian coordinate advection
	I'll still be responsible for setting some beta^i_,t gauge
	so for the L_beta Lie derivative, we have some options:
	1) none (for no-shift)
	2) finite difference
	3) finite volume / absorb into the eigensystem
	4) Lagrangian coordinates
	and this should be a separate variable, separate of the shift gauge

	so 
	one variable for what beta^i_,t is
	another variable for how to 
	--]=]
--]]
function Z4_2004Bona:init(args)
	local solver = args.solver

	local fluxVars = table{
		{name='a_l', type='real3'},
		{name='d_lll', type='_3sym3'},
		{name='K_ll', type='sym3'},
		{name='Theta', type='real'},
		{name='Z_l', type='real3'},
	}

	self.consVars = table{
		{name='alpha', type='real'},
		{name='gamma_ll', type='sym3'},
	}:append(fluxVars)


	--[[
	how are shift conditions impmlemented?
	options for determining beta^i:
	1) by solving a constraint equation (minimal distortion elliptic solves it with a numerical Poisson solver)
	2) by solving an initial value problem

	Once beta^i is determined, how are the variables iterated?
	This question is split into a) and b):
	a) How are the source-only variables iterated wrt beta^i?
	options:
	1) put the Lie derivative terms into the source side of these variables 
	2) give them -beta^i eigenvalues?  this is the equivalent of rewriting the hyperbolic vars associated with these (a_i,d_kij) back into first-derivative 0th order vars (alpha, gamma_ij)

	b) How are the flux variables iterated wrt beta^i?
	options:
	1) this would be solved by offsetting the eigenvalues 
		with the offset eigenvalues, 
		the hyperbolic state vars' contributions get incorporated into the flux,
		but there are still some shift-based terms that end up in the source ...
		... so the shift is split between the flux and source ...
	2) ... and adding a few terms to the source
	--]]

	self.useShift = args.useShift or 'none'
	
	-- set to false to disable rho, S_i, S^ij
	self.useStressEnergyTerms = true

	if self.useShift ~= 'none' then
		self.consVars:insert{name='beta_u', type='real3'}

		if self.useShift == 'MinimalDistortionElliptic' 
		or self.useShift == 'MinimalDistortionEllipticEvolve' 
		then
			self.consVars:insert{name='betaLap_u', type='real3'}
		end
	end


	--[[
	solve a smaller eigendecomposition that doesn't include the rows of variables whose d/dt is zero.
	kind of like how ideal MHD is an 8-var system, but the flux jacobian solved is a 7x7 because Bx,t = 0.
	TODO make this a ctor arg - so solvers can run in parallel with and without this
	...or why don't I just scrap the old code, because this runs a lot faster.
	--]]
	self.noZeroRowsInFlux = true
	--self.noZeroRowsInFlux = false
	if args.noZeroRowsInFlux ~= nil then
		self.noZeroRowsInFlux = args.noZeroRowsInFlux
	end

	-- NOTE this doesn't work when using shift ... because then all the eigenvalues are -beta^i, so none of them are zero (except the source-only alpha, beta^i, gamma_ij)
	-- with the exception of the lagrangian shift.  if we split that operator out then we can first solve the non-shifted system, then second advect it by the shift vector ...
	--if self.useShift ~= 'none' then
	--	self.noZeroRowsInFlux = false
	--end

	if not self.noZeroRowsInFlux then
		-- skip alpha and gamma
		self.numWaves = Struct.countScalars{vars=fluxVars}
		assert(self.numWaves == 30)
	else
		-- ok Z4 has a_x, d_xij, K_ij, Theta, Z_i ...
		-- which is 17 waves
		-- the min/max which HLL uses are the same waves iirc
		-- but... I should fix that. 
		self.numWaves = 13
	end

	-- only count int vars after the shifts have been added
	self:cdefAllVarTypes(solver, self.consVars)	-- have to call before countScalars in eqn:init
	self.numIntStates = Struct.countScalars{vars=self.consVars}

	-- now add in the source terms (if you want them)
	if self.useStressEnergyTerms then
		self.consVars:append{
			--stress-energy variables:
			{name='rho', type='real'},					--1: n_a n_b T^ab
			{name='S_u', type='real3'},				--3: -gamma^ij n_a T_aj
			{name='S_ll', type='sym3'},				--6: gamma_i^c gamma_j^d T_cd
		}								
	end
	self.consVars:append{
		--constraints:              
		{name='H', type='real'},					--1
		{name='M_u', type='real3'},				--3
	}

	self.eigenVars = table{
		{name='alpha', type='real'},
		{name='alpha_sqrt_f', type='real'},
		{name='gamma_ll', type='sym3'},
		{name='gamma_uu', type='sym3'},
		-- sqrt(gamma^jj) needs to be cached, otherwise the Intel kernel stalls (for seconds on end)
		{name='sqrt_gammaUjj', type='real3'},
	}

	-- hmm, only certain shift methods actually use beta_u ...
	if self.useShift ~= 'none' then
		self.eigenVars:insert{name='beta_u', type='real3'}
	end


	-- build stuff around consVars	
	Z4_2004Bona.super.init(self, args)


	if self.useShift == 'MinimalDistortionElliptic' then
		local MinimalDistortionEllipticShift = require 'hydro.op.gr-shift-mde'
		self.solver.ops:insert(MinimalDistortionEllipticShift{solver=self.solver})
	elseif self.useShift == 'LagrangianCoordinates' then
		local LagrangianCoordinateShift = require 'hydro.op.gr-shift-lc'
		self.solver.ops:insert(LagrangianCoordinateShift{solver=self.solver})
	end
end

function Z4_2004Bona:createInitState()
	Z4_2004Bona.super.createInitState(self)
	self:addGuiVars{
		
		-- from 2004 Bona et al, "A symmetry breaking..." eqn A.20
		{name='m', value=2},
	
		-- convergence between finite-difference of alpha,i and alpha a_i
		{name='a_convCoeff', value=0},
		
		-- convergence between finite-difference of 1/2 gamma_ij,k and d_kij
		{name='d_convCoeff', value=0},
	
	}
	-- TODO add shift option
	-- but that means moving the consVars construction to the :init()
end

function Z4_2004Bona:initCodeModules()
	Z4_2004Bona.super.initCodeModules(self)
	local solver = self.solver
	
	solver.modules:add{
		name = 'calc_gamma_ll',
		code = [[
#define calc_gamma_ll(U, x)	((U)->gamma_ll)
]],
	}

	solver.modules:add{
		name = 'calc_gamma_uu',
		depends = {
			'eqn.cons_t',
		},
		code = self:template[[
sym3 calc_gamma_uu(const global <?=eqn.cons_t?>* U, real3 x) {
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	return gamma_uu;
}
]],
	}
end

function Z4_2004Bona:initCodeModule_calcDT()
	local solver = self.solver
	solver.modules:add{
		name = 'calcDT',
		depends = {
			'eqn.cons_t',
			'initCond.codeprefix',	-- calc_f
		},
		code = self:template[[
kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
		
	const global <?=eqn.cons_t?>* U = UBuf + index;
	
	//the only advantage of this calcDT over the default is that here this sqrt(f) and det(gamma_ij) is only called once
	real f_alphaSq = calc_f_alphaSq(U->alpha);
	real det_gamma = sym3_det(U->gamma_ll);
	real alpha_sqrt_f = sqrt(f_alphaSq);

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		
		<? if side == 0 then ?>
		real gammaUjj = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;
		<? elseif side == 1 then ?>
		real gammaUjj = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;
		<? elseif side == 2 then ?>
		real gammaUjj = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;
		<? end ?>	
		real lambdaLight = U->alpha * sqrt(gammaUjj);
		
		real lambdaGauge = sqrt(gammaUjj)* alpha_sqrt_f;
		real lambda = (real)max(lambdaGauge, lambdaLight);

		<? if eqn.useShift ~= 'none' then ?>
		real betaUi = U->beta_u.s<?=side?>;
		<? else ?>
		const real betaUi = 0.;
		<? end ?>
		
		real lambdaMin = (real)min((real)0., -betaUi - lambda);
		real lambdaMax = (real)max((real)0., -betaUi + lambda);
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
		absLambdaMax = max((real)1e-9, absLambdaMax);
		dt = (real)min(dt, solver->grid_dx.s<?=side?> / absLambdaMax);
	}<? end ?>
	dtBuf[index] = dt; 
}
]],
	}
end

function Z4_2004Bona:initCodeModule_fluxFromCons()
	self.solver.modules:add{
		name = 'fluxFromCons',
		depends = {
			'eqn.cons_t',
			'solver.solver_t',
			'normal_t',
		},
		code = self:template[[

<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	real f_alpha = calc_f_alpha(U.alpha);
	
	real det_gamma = sym3_det(U.gamma_ll);
	sym3 gamma_uu = sym3_inv(U.gamma_ll, det_gamma);
	
	real alpha = U.alpha;
	real Theta = U.Theta;
	
	real3 Z_l = U.Z_l;
	real3 a_l = U.a_l;
	sym3 gamma_ll = U.gamma_ll;
	sym3 K_ll = U.K_ll;
	_3sym3 d_lll = U.d_lll;

	if (false) {}
	<? for side=0,solver.dim-1 do ?>
	else if (n.side == <?=side?>) {
		Z_l = real3_swap<?=side?>(Z_l);
		a_l = real3_swap<?=side?>(a_l);
		gamma_ll = sym3_swap<?=side?>(gamma_ll);
		K_ll = sym3_swap<?=side?>(K_ll);
		d_lll = _3sym3_swap<?=side?>(d_lll);
		gamma_uu = sym3_swap<?=side?>(gamma_uu);
	}
	<? end ?>
	else {
		alpha = 0./0.;
		Theta = 0./0.;
<? for k,xk in ipairs(xNames) do
?>		a_l.<?=xk?> = 0./0.;
		Z_l.<?=xk?> = 0./0.;
<? end
?>	
<? for ij,xij in ipairs(symNames) do
?>		K_ll.<?=xij?> = 0./0.;
		gamma_ll.<?=xij?> = 0./0.;
<? end
?>
<? for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>		d_lll.<?=xk?>.<?=xij?> = 0./0.;
<?	end
end
?>	}

	<?=eqn.cons_t?> F = {.ptr={ 0./0. }};
	//for (int i = 0; i < numStates; ++i) {
	//	F.ptr[i] = 0.;
	//}


	// BEGIN CUT from numerical-relativity-codegen/flux_matrix_output/z4_noZeroRows.html
	real tmp1 = K_ll.xy * gamma_uu.xy;
	real tmp2 = K_ll.xz * gamma_uu.xz;
	real tmp3 = K_ll.yz * gamma_uu.yz;
	real tmp4 = 2. * tmp3;
	real tmp6 = K_ll.yy * gamma_uu.yy;
	real tmp7 = K_ll.zz * gamma_uu.zz;
	real tmp103 = gamma_uu.xx * gamma_uu.yy;
	real tmp106 = gamma_uu.xy * gamma_uu.xy;
	real tmp107 = gamma_uu.xx * gamma_uu.yz;
	real tmp111 = gamma_uu.xy * gamma_uu.xz;
	real tmp113 = gamma_uu.xx * gamma_uu.zz;
	real tmp116 = gamma_uu.xz * gamma_uu.xz;
	real tmp125 = gamma_uu.xy * gamma_uu.yz;
	real tmp128 = gamma_uu.xz * gamma_uu.yy;
	real tmp129 = gamma_uu.xy * gamma_uu.zz;
	real tmp132 = gamma_uu.xz * gamma_uu.yz;
	F.alpha = 0.;
	F.gamma_ll.xx = 0.;
	F.gamma_ll.xy = 0.;
	F.gamma_ll.xz = 0.;
	F.gamma_ll.yy = 0.;
	F.gamma_ll.yz = 0.;
	F.gamma_ll.zz = 0.;
	F.a_l.x = f_alpha * (tmp6 + tmp7 + K_ll.xx * gamma_uu.xx + tmp4 - 2. * Theta + 2. * tmp2 + 2. * tmp1);
	F.a_l.y = 0.;
	F.a_l.z = 0.;
	F.d_lll.x.xx = K_ll.xx * alpha;
	F.d_lll.x.xy = K_ll.xy * alpha;
	F.d_lll.x.xz = K_ll.xz * alpha;
	F.d_lll.x.yy = K_ll.yy * alpha;
	F.d_lll.x.yz = K_ll.yz * alpha;
	F.d_lll.x.zz = K_ll.zz * alpha;
	F.d_lll.y.xx = 0.;
	F.d_lll.y.xy = 0.;
	F.d_lll.y.xz = 0.;
	F.d_lll.y.yy = 0.;
	F.d_lll.y.yz = 0.;
	F.d_lll.y.zz = 0.;
	F.d_lll.z.xx = 0.;
	F.d_lll.z.xy = 0.;
	F.d_lll.z.xz = 0.;
	F.d_lll.z.yy = 0.;
	F.d_lll.z.yz = 0.;
	F.d_lll.z.zz = 0.;
	F.K_ll.xx = -alpha * (2. * d_lll.z.xz * gamma_uu.zz + 2. * d_lll.z.xy * gamma_uu.yz + d_lll.z.xx * gamma_uu.xz + 2. * d_lll.y.xz * gamma_uu.yz + 2. * d_lll.y.xy * gamma_uu.yy + d_lll.y.xx * gamma_uu.xy + 2. * Z_l.x - a_l.x - d_lll.x.yy * gamma_uu.yy - 2. * d_lll.x.yz * gamma_uu.yz - d_lll.x.zz * gamma_uu.zz);
	F.K_ll.xy = (-alpha * (2. * d_lll.z.yz * gamma_uu.zz + 2. * d_lll.z.yy * gamma_uu.yz + 2. * d_lll.y.yz * gamma_uu.yz + 2. * d_lll.y.yy * gamma_uu.yy + 2. * d_lll.x.yz * gamma_uu.xz + 2. * d_lll.x.yy * gamma_uu.xy + 2. * Z_l.y - a_l.y)) / 2.;
	F.K_ll.xz = (-alpha * (2. * d_lll.z.zz * gamma_uu.zz + 2. * d_lll.z.yz * gamma_uu.yz + 2. * d_lll.y.zz * gamma_uu.yz + 2. * d_lll.y.yz * gamma_uu.yy + 2. * d_lll.x.zz * gamma_uu.xz + 2. * d_lll.x.yz * gamma_uu.xy + 2. * Z_l.z - a_l.z)) / 2.;
	F.K_ll.yy = alpha * (d_lll.z.yy * gamma_uu.xz + d_lll.y.yy * gamma_uu.xy + d_lll.x.yy * gamma_uu.xx);
	F.K_ll.yz = alpha * (d_lll.z.yz * gamma_uu.xz + d_lll.y.yz * gamma_uu.xy + d_lll.x.yz * gamma_uu.xx);
	F.K_ll.zz = alpha * (d_lll.z.zz * gamma_uu.xz + d_lll.y.zz * gamma_uu.xy + d_lll.x.zz * gamma_uu.xx);
	F.Theta = -alpha * (d_lll.z.yz * tmp129 - d_lll.z.yz * tmp132 + d_lll.z.yy * tmp125 - d_lll.z.yy * tmp128 + d_lll.z.xz * tmp113 - d_lll.z.xz * tmp116 + d_lll.z.xy * tmp107 - d_lll.z.xy * tmp111 + d_lll.y.zz * tmp132 + d_lll.y.yz * tmp128 - d_lll.y.zz * tmp129 + d_lll.y.xz * tmp107 - d_lll.y.xz * tmp111 - d_lll.y.yz * tmp125 + d_lll.y.xy * tmp103 - d_lll.y.xy * tmp106 + d_lll.x.zz * tmp116 + 2. * d_lll.x.yz * tmp111 - d_lll.x.zz * tmp113 + d_lll.x.yy * tmp106 - 2. * d_lll.x.yz * tmp107 + Z_l.z * gamma_uu.xz - d_lll.x.yy * tmp103 + Z_l.y * gamma_uu.xy + Z_l.x * gamma_uu.xx);
	F.Z_l.x = alpha * (tmp1 + tmp2 + tmp6 + tmp4 + tmp7 - Theta);
	F.Z_l.y = -alpha * (K_ll.yz * gamma_uu.xz + K_ll.yy * gamma_uu.xy + K_ll.xy * gamma_uu.xx);
	F.Z_l.z = -alpha * (K_ll.zz * gamma_uu.xz + K_ll.yz * gamma_uu.xy + K_ll.xz * gamma_uu.xx);	
	// END CUT

	if (false) {}
	<? for side=0,solver.dim-1 do ?>
	else if (n.side == <?=side?>) {
		F.Z_l = real3_swap<?=side?>(F.Z_l);
		F.a_l = real3_swap<?=side?>(F.a_l);
		F.gamma_ll = sym3_swap<?=side?>(F.gamma_ll);
		F.K_ll = sym3_swap<?=side?>(F.K_ll);
		F.d_lll = _3sym3_swap<?=side?>(F.d_lll);
	}
	<? end ?>
	else {
		alpha = 0./0.;
		Theta = 0./0.;
<? for k,xk in ipairs(xNames) do
?>		a_l.<?=xk?> = 0./0.;
		Z_l.<?=xk?> = 0./0.;
<? end
?>	
<? for ij,xij in ipairs(symNames) do
?>		K_ll.<?=xij?> = 0./0.;
		gamma_ll.<?=xij?> = 0./0.;
<? end
?>
<? for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>		d_lll.<?=xk?>.<?=xij?> = 0./0.;
<?	end
end
?>	}


<? 
if eqn.useShift ~= 'none' then
?>	//beta^i_,t = 0 + source terms
	F.beta_u = real3_zero;
<?	if self.useShift == 'MinimalDistortionElliptic' 
	or self.useShift == 'MinimalDistortionEllipticEvolve' 
	then
?>	F.betaLap_u = real3_zero;
<?	end
end 
?>

	return F;
}
]],
	}
end

function Z4_2004Bona:initCodeModule_setFlatSpace()
	self.solver.modules:add{
		name = 'setFlatSpace',
		depends = {
			'solver.solver_t',
			'eqn.cons_t',
		},
		code = self:template[[
void setFlatSpace(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* U,
	real3 x
) {
	U->alpha = 1.;
	U->gamma_ll = sym3_ident;
	U->a_l = real3_zero;
	U->d_lll.x = sym3_zero;
	U->d_lll.y = sym3_zero;
	U->d_lll.z = sym3_zero;
	U->K_ll = sym3_zero;
	U->Theta = 0.;
	U->Z_l = real3_zero;
<? if eqn.useShift ~= 'none' then 
?>	U->beta_u = real3_zero;
<? end 
?>	
<? if eqn.useStressEnergyTerms then ?>
	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	U->H = 0;
	U->M_u = real3_zero;
}
]],
	}
end

function Z4_2004Bona:getModuleDependsSolver()
	return table(Z4_2004Bona.super.getModuleDependsSolver(self))
	:append{
		'rotate',	-- real3_swap ... though sym3_swap and _3sym3_swap are in their respective modules ...
		'initCond.codeprefix',	-- calc_f
	}
end

function Z4_2004Bona:getModuleDependsApplyInitCond()
	return {
		'coordMap',
		'coord_g_ll',
		'rescaleFromCoord/rescaleToCoord',
	}
end


-- always true, except for initAnalytical, but I don't have that coded yet
Z4_2004Bona.needsInitDerivs = true

function Z4_2004Bona:getInitCondCode()
	if self.initCond.initAnalytical then
		error("TODO - can't handle analytical initial conditions yet")
	end
	
	if self.initCond.useBSSNVars then
		return self:template([[
kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cellBuf[index].pos;
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real alpha = 1.;
	real W = 1.;
	real K = 0.;
	real3 LambdaBar_U = real3_zero;
	real3 beta_U = real3_zero;
	real3 B_U = real3_zero;
	sym3 epsilon_LL = sym3_zero;
	sym3 ABar_LL = sym3_zero;

	real rho = 0.;

	<?=code?>

	//for (int i = 0; i < numStates; ++i) {
	//	U->ptr[i] = 0.;
	//}
	*U = (<?=eqn.cons_t?>){.ptr={0}};

	U->alpha = alpha;

	// gammaHat_IJ = delta_IJ
	// gamma_ij = e_i^I e_j^J (epsilon_IJ + gammaHat_IJ) / W^2
	sym3 gammaBar_LL = sym3_add(epsilon_LL, sym3_ident);
	sym3 gamma_LL = sym3_real_mul(gammaBar_LL, 1. / (W*W));
	U->gamma_ll = sym3_rescaleToCoord_LL(gamma_LL, x);
	
	// K_ij = e_i^I e_j^J (ABar_IJ + gammaBar_IJ K/3) / W^2
	U->K_ll = sym3_rescaleToCoord_LL(
		sym3_add(
			sym3_real_mul(ABar_LL, 1. / (W*W)),
			sym3_real_mul(gamma_LL, K / 3.)
		), x);

	U->Theta = 0.;
	U->Z_l = real3_zero;

<? if eqn.useShift ~= 'none' then
?>	U->beta_u = real3_rescaleFromCoord_U(beta_U);
<? end -- TODO support for hyperbolic gamma driver, so we can read B_U
?>

<? if eqn.useStressEnergyTerms then ?>
	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	U->H = 0;
	U->M_u = real3_zero;
}

kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	global <?=solver.coord.cell_t?> const *cellBuf 
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

<? 
for i=1,solver.dim do 
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = (
		log(U[solver->stepsize.<?=xi?>].alpha) 
		- log(U[-solver->stepsize.<?=xi?>].alpha)
	) / (2. * solver->grid_dx.s<?=i-1?>);
	<? for jk,xjk in ipairs(symNames) do ?>
	U->d_lll.<?=xi?>.<?=xjk?> = .5 * (
		U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> 
		- U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>
	) / (2. * solver->grid_dx.s<?=i-1?>);
	<? end ?>
<? 
end 
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = 0;
	U->d_lll.<?=xi?> = sym3_zero;
<?
end
?>
}
]], 	{
			code = self.initCond:getInitCondCode(self.solver),
		})

	end

	return self:template([[
kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = coord_g_ll(x);
	sym3 K_ll = sym3_zero;

	//TODO more stress-energy vars 
	real rho = 0.;

	<?=code?>

	U->alpha = alpha;
	U->gamma_ll = gamma_ll;
	U->K_ll = K_ll;
	U->Theta = 0.;
	U->Z_l = real3_zero;

<? if eqn.useShift ~= 'none' then
?>	U->beta_u = beta_u;
<? end
?>

<? if eqn.useStressEnergyTerms then ?>
	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
<? end ?>
	U->H = 0;
	U->M_u = real3_zero;
}

kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	global <?=solver.coord.cell_t?> const *cellBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);

<? 
for i=1,solver.dim do 
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = (
		log(U[solver->stepsize.<?=xi?>].alpha) 
		- log(U[-solver->stepsize.<?=xi?>].alpha)
	) / (2. * solver->grid_dx.s<?=i-1?>);
	<? for jk,xjk in ipairs(symNames) do ?>
	U->d_lll.<?=xi?>.<?=xjk?> = .5 * (
		U[solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?> 
		- U[-solver->stepsize.<?=xi?>].gamma_ll.<?=xjk?>
	) / (2. * solver->grid_dx.s<?=i-1?>);
	<? end ?>
<? 
end 
for i=solver.dim+1,3 do
	local xi = xNames[i]
?>
	U->a_l.<?=xi?> = 0;
	U->d_lll.<?=xi?> = sym3_zero;
<?
end
?>
}
]], {
		code = self.initCond:getInitCondCode(self.solver),
	})
end

Z4_2004Bona.solverCodeFile = 'hydro/eqn/z4.cl'

Z4_2004Bona.predefinedDisplayVars = {
	'U alpha',
--[[ for x dir only
	'U gamma_ll x x',
	'U d_lll_x x x',
	'U K_ll x x',
	'U Theta',
	'U Z_l x',
--]]	
-- [[ for all dirs
	'U gamma_ll norm',
	'U a_l mag',
	-- TODO a 3sym3 Frobenius norm - use for 'd'?
	'U d_lll x norm',
	'U d_lll y norm',
	'U d_lll z norm',
	'U K_ll norm',
	'U Theta',
	'U Z_l mag',
--]]
	'U H',
	'U M_u mag',
	'U volume',
	'U f',
}

function Z4_2004Bona:getDisplayVars()
	local vars = Z4_2004Bona.super.getDisplayVars(self)

	vars:append{
		{name='volume', code='value.vreal = U->alpha * sqrt(sym3_det(U->gamma_ll));'},
		{name='f', code='value.vreal = calc_f(U->alpha);'},
		{name='df/dalpha', code='value.vreal = calc_dalpha_f(U->alpha);'},
		
		-- is expansion really just -K?
		{name='expansion', code=[[
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	value.vreal = -sym3_dot(gamma_uu, U->K_ll);
]]		},
	}:append{
--[=[
--[[
Alcubierre 3.1.1
Baumgarte & Shapiro 2.90
H = R + K^2 - K^ij K_ij - 16 pi rho
for rho = n_a n_b T^ab (B&S eqn 2.89)
and n_a = -alpha t_,a (B&S eqns 2.19, 2.22, 2.24)

momentum constraints
--]]
		{H = [[
	.5 * 
]]		},
--]=]
	}

	-- shift-less gravity only
	-- gravity with shift is much more complex
	-- TODO add shift influence (which is lengthy)
	vars:insert{name='gravity', code=[[
	real det_gamma = sym3_det(U->gamma_ll);
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	value.vreal3 = real3_real_mul(sym3_real3_mul(gamma_uu, U->a_l), -U->alpha * U->alpha);
]], type='real3'}

	vars:insert{name='alpha vs a_i', code=template([[
	if (OOB(1,1)) {
		value.vreal3 = real3_zero;
	} else {
		<? for i=1,solver.dim do
			local xi = xNames[i]
		?>{
			real di_log_alpha = (
				log(U[solver->stepsize.<?=xi?>].alpha)
				- log(U[-solver->stepsize.<?=xi?>].alpha)
			) / (2. * solver->grid_dx.s<?=i-1?>);
			value.vreal3.<?=xi?> = fabs(di_log_alpha - U->a_l.<?=xi?>);
		}<? end ?>
		<? for i=solver.dim+1,3 do
			local xi = xNames[i]
		?>{
			value.vreal3.<?=xi?> = 0;
		}<? end ?>
	}
]], {
	solver = self.solver,
	xNames = xNames,
}), type='real3'}

	-- d_kij = gamma_ij,k
	for i,xi in ipairs(xNames) do
		vars:insert{name='gamma_ij vs d_'..xi..'ij', code=template([[
	if (OOB(1,1)) {
		value.vsym3 = sym3_zero;
	} else {
		<? if i <= solver.dim then ?>
		sym3 di_gamma_jk = sym3_real_mul(
			sym3_sub(
				U[solver->stepsize.<?=xi?>].gamma_ll, 
				U[-solver->stepsize.<?=xi?>].gamma_ll
			), 
			.5 / (2. * solver->grid_dx.s<?=i-1?>)
		);
		<? else ?>
		sym3 di_gamma_jk = sym3_zero;
		<? end ?>
		value.vsym3 = sym3_sub(di_gamma_jk, sym3_real_mul(U->d_lll.<?=xi?>, 2.));
		value.vsym3 = (sym3){<?
	for jk,xjk in ipairs(symNames) do 
?>			.<?=xjk?> = fabs(value.vsym3.<?=xjk?>),
<?	end
?>		};
	}
]], {
	i = i,
	xi = xi,
	xNames = xNames,
	symNames = symNames,
	solver = self.solver,
}), type='sym3'}
	end

	return vars
end

function Z4_2004Bona:eigenWaveCodePrefix(n, eig, x, waveIndex)
	return template([[
	real eig_lambdaLight = <?=eig?>.sqrt_gammaUjj.s[<?=n?>.side] * <?=eig?>.alpha;
	real eig_lambdaGauge = <?=eig?>.sqrt_gammaUjj.s[<?=n?>.side] * <?=eig?>.alpha_sqrt_f;
]], {
		eig = '('..eig..')',
		n = n,
	})
end

function Z4_2004Bona:eigenWaveCode(n, eig, x, waveIndex)
	-- TODO find out if -- if we use the lagrangian coordinate shift operation -- do we still need to offset the eigenvalues by -beta^i?
	--local shiftingLambdas = self.useShift ~= 'none'
	--and self.useShift ~= 'LagrangianCoordinates'

	local betaUi
	if self.useShift ~= 'none' then
		betaUi = '('..eig..').beta_u.s['..n..'.side]'
	else
		betaUi = '0'
	end

	if not self.noZeroRowsInFlux then
		if waveIndex == 0 then
			return '-'..betaUi..' - eig_lambdaGauge'
		elseif waveIndex >= 1 and waveIndex <= 6 then
			return '-'..betaUi..' - eig_lambdaLight'
		elseif waveIndex >= 7 and waveIndex <= 23 then
			return '-'..betaUi
		elseif waveIndex >= 24 and waveIndex <= 29 then
			return '-'..betaUi..' + eig_lambdaLight'
		elseif waveIndex == 30 then
			return '-'..betaUi..' + eig_lambdaGauge'
		end
	else	-- noZeroRowsInFlux 
		-- noZeroRowsInFlux implies useShift == 'none'
		if waveIndex == 0 then
			return '-'..betaUi..' - eig_lambdaGauge'
		elseif waveIndex >= 1 and waveIndex <= 5 then
			return '-'..betaUi..' - eig_lambdaLight'
		elseif waveIndex == 6 then
			return '-'..betaUi
		elseif waveIndex >= 7 and waveIndex <= 11 then
			return '-'..betaUi..' + eig_lambdaLight'
		elseif waveIndex == 12 then
			return '-'..betaUi..' + eig_lambdaGauge'
		end
	end
	error'got a bad waveIndex'
end

function Z4_2004Bona:eigenWaveMinCode(n, eig, x)
	local betaUi
	if self.useShift ~= 'none' then
		betaUi = '('..eig..').beta_u.s['..n..'.side]'
	else
		betaUi = '0'
	end
	return 'min(-'..betaUi..', min(-'..betaUi..' - eig_lambdaGauge, -'..betaUi..' - eig_lambdaLight))'
end

function Z4_2004Bona:eigenWaveMaxCode(n, eig, x)
	local betaUi
	if self.useShift ~= 'none' then
		betaUi = '('..eig..').beta_u.s['..n..'.side]'
	else
		betaUi = '0'
	end
	return 'max(-'..betaUi..', max(-'..betaUi..' + eig_lambdaGauge, -'..betaUi..' + eig_lambdaLight))'
end

function Z4_2004Bona:consWaveCodePrefix(n, U, x, waveIndex)
	return template([[
	real det_gamma = sym3_det(<?=U?>.gamma_ll);
	sym3 gamma_uu = sym3_inv(<?=U?>.gamma_ll, det_gamma);
	real sqrt_gammaUjj = 0./0.;
	if (<?=n?>.side == 0) {
		sqrt_gammaUjj = sqrt(gamma_uu.xx);
	} else if (<?=n?>.side == 1) {
		sqrt_gammaUjj = sqrt(gamma_uu.yy);
	} else if (<?=n?>.side == 2) {
		sqrt_gammaUjj = sqrt(gamma_uu.zz);
	}
	real eig_lambdaLight = <?=U?>.alpha * sqrt_gammaUjj;
	real f_alphaSq = calc_f_alphaSq(<?=U?>.alpha);
	real eig_lambdaGauge = sqrt(f_alphaSq) * sqrt_gammaUjj;
]], {
		U = '('..U..')',
		n = n,
	})
end
Z4_2004Bona.consWaveCode = Z4_2004Bona.eigenWaveCode

return Z4_2004Bona

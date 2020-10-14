--[[
based on whatever my numerical_relativity_codegen z4 is based on, which is probably a Bona-Masso paper,
probably 2004 Bona et al "A symmetry-breaking mechanism for the Z4 general-covariant evolution system"
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
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
		-- sqrt(n_i n_j gamma^ij) needs to be cached, otherwise the Intel kernel stalls (for seconds on end)
		{name='sqrt_gammaUnn', type='real'},
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

	-- used by hll, roe, weno, plm ... anything that uses eigenvalues or eigenvector transforms
	solver.modules:add{
		name = 'eigen_forInterface',
		depends = {
			'solver.solver_t',
			'eqn.eigen_t',
			'eqn.cons_t',
			'normal_t',
		},
		code = self:template[[
//used for interface eigen basis
<?=eqn.eigen_t?> eigen_forInterface(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> UL,
	<?=eqn.cons_t?> UR,
	real3 x,
	normal_t n
) {
	real alpha = .5 * (UL.alpha + UR.alpha);
	sym3 avg_gamma = (sym3){
		.xx = .5 * (UL.gamma_ll.xx + UR.gamma_ll.xx),
		.xy = .5 * (UL.gamma_ll.xy + UR.gamma_ll.xy),
		.xz = .5 * (UL.gamma_ll.xz + UR.gamma_ll.xz),
		.yy = .5 * (UL.gamma_ll.yy + UR.gamma_ll.yy),
		.yz = .5 * (UL.gamma_ll.yz + UR.gamma_ll.yz),
		.zz = .5 * (UL.gamma_ll.zz + UR.gamma_ll.zz),
	};
	real det_avg_gamma = sym3_det(avg_gamma);

	<?=eqn.eigen_t?> eig;
	eig.alpha = alpha;
	eig.alpha_sqrt_f = sqrt(calc_f_alphaSq(alpha));
	eig.gamma_uu = sym3_inv(avg_gamma, det_avg_gamma);

<? if solver.coord.vectorComponent == 'cartesian' then ?>
// I'm using .side for holonomic(coordinate) and anholonomic(orthonormal)
//but for cartesian vector componets there is no .side, just .n, which is covariant iirc
//and I haven't derived the flux in arbitrary-normal form, just in x-axis form (and i swap x<->y or z to calculate their fluxes)
//so here I'm going to just wing it
	real3 n_l = normal_l1(n);
	real gammaUnn = real3_weightedLenSq(n_l, eig.gamma_uu);
<? else ?>
	real gammaUnn = 0./0.;
	if (n.side == 0) {
		gammaUnn = eig.gamma_uu.xx;
	} else if (n.side == 1) {
		gammaUnn = eig.gamma_uu.yy;
	} else if (n.side == 2) {
		gammaUnn = eig.gamma_uu.zz;
	}
<? end ?>

	eig.sqrt_gammaUnn = sqrt(gammaUnn);
	
<? if eqn.useShift ~= 'none' then ?>
	eig.beta_u = real3_real_mul(real3_add(UL.beta_u, UR.beta_u), .5);
<? end ?>

	return eig;
}
]],
	}

	-- not used anymore, replaced in calcDT by eqn:consMinWaveCode/eqn:consMaxWaveCode eigenvalue inlining
	solver.modules:add{
		name = 'calcCellMinMaxEigenvalues',
		depends = {
			'range_t',
			'normal_t',
			'eqn.cons_t',
		},
		code = self:template[[
range_t calcCellMinMaxEigenvalues(
	const global <?=eqn.cons_t?>* U,
	real3 x,
	normal_t n
) {
	real det_gamma = sym3_det(U->gamma_ll);

<? if solver.coord.vectorComponent == 'cartesian' then ?>
	sym3 gamma_uu = sym3_inv(U->gamma_ll, det_gamma);
	real3 n_l = normal_l1(n);
	real gammaUnn = real3_weightedLenSq(n_l, gamma_uu);
<? else ?>
	real gammaUnn = 0./0.;
	if (n.side == 0) {
		gammaUnn = (U->gamma_ll.yy * U->gamma_ll.zz - U->gamma_ll.yz * U->gamma_ll.yz) / det_gamma;
	} else if (n.side == 1) {
		gammaUnn = (U->gamma_ll.xx * U->gamma_ll.zz - U->gamma_ll.xz * U->gamma_ll.xz) / det_gamma;
	} else if (n.side == 2) {
		gammaUnn = (U->gamma_ll.xx * U->gamma_ll.yy - U->gamma_ll.xy * U->gamma_ll.xy) / det_gamma;
	}
<? end ?>
	real sqrt_gammaUnn = sqrt(gammaUnn);
	real lambdaLight = U->alpha * sqrt_gammaUnn;
	
	real f_alphaSq = calc_f_alphaSq(U->alpha);
	real lambdaGauge = sqrt(f_alphaSq) * sqrt_gammaUnn;

	real lambdaMax = max(lambdaGauge, lambdaLight);
	real lambdaMin = -lambdaMin;

	<? if eqn.useShift ~= 'none' then ?>
	lambdaMin -= normal_vecDotN1(n, U->beta_u);
	lambdaMax -= normal_vecDotN1(n, U->beta_u);
	<? end ?>

	return (range_t){
		.min = lambdaMin, 
		.max = lambdaMax,
	};
}
]],
	}

	-- used by plm
	solver.modules:add{
		name = 'eigen_forCell',
		depends = {
			'solver.solver_t',
			'eqn.cons_t',
			'normal_t',
			'initCond.codeprefix',	-- calc_f
		},
		code = self:template[[
//used by PLM, and by the default fluxFromCons (used by hll, or roe when roeUseFluxFromCons is set)
<?=eqn.eigen_t?> eigen_forCell(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	<?=eqn.eigen_t?> eig;
	eig.alpha = U.alpha;
	eig.alpha_sqrt_f = sqrt(calc_f_alphaSq(U.alpha));
	real det_gamma = sym3_det(U.gamma_ll);
	eig.gamma_uu = sym3_inv(U.gamma_ll, det_gamma);

<? if solver.coord.vectorComponent == 'cartesian' then ?>
	real3 n_l = normal_l1(n);
	real gammaUnn = real3_weightedLenSq(n_l, eig.gamma_uu);
<? else ?>
	real gammaUnn = 0./0.;
	if (n.side == 0) {
		gammaUnn = eig.gamma_uu.xx;
	} else if (n.side == 1) {
		gammaUnn = eig.gamma_uu.yy;
	} else if (n.side == 2) {
		gammaUnn = eig.gamma_uu.zz;
	}
<? end ?>

	eig.sqrt_gammaUnn = sqrt(gammaUnn);
	
	<? if eqn.useShift ~= 'none' then ?>
	eig.beta_u = U.beta_u;
	<? end ?>

	return eig;
}
]],
	}

	-- used by roe, weno, some plm
	solver.modules:add{
		name = 'eigen_left/rightTransform',
		depends = {
			'solver.solver_t',
			'eqn.cons_t',
			'normal_t',
		},
		code = self:template[[
//TODO these were based no noZeroRowsInFlux==false (I think) so maybe/certainly they are out of date
<?=eqn.waves_t?> eigen_leftTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> inputU,
	real3 x,
	normal_t n
) {
	<?=eqn.waves_t?> results = (<?=eqn.waves_t?>){.ptr={0. / 0.}};
	//for (int i = 0; i < numWaves; ++i) {
	//	results.ptr[i] = 0;
	//}

#error here
#if 0	//don't enable this.  it's made for > waves than I'm using, so it will cause buffer corruption

	real3 a_l = real3_swap(inputU.a_l, n.side);							//0-2
	_3sym3 d_lll = _3sym3_swap(inputU.d_lll, n.side);					//3-20 ... .x = 3-8, .y = 9-14, .z = 15-20
	sym3 K_ll = sym3_swap(inputU.K_ll, n.side);							//21-26
	real Theta = inputU.Theta;												//27
	real3 Z_l = real3_swap(inputU.Z_l, n.side);							//28-30
	
	sym3 gamma_ll = sym3_swap(eig.gamma_ll, n.side);
	sym3 gamma_uu = sym3_swap(eig.gamma_uu, n.side);

	real sqrt_gammaUUxx = sqrt(gamma_uu.xx);
	real gammaUUxx_toThe_3_2 = sqrt_gammaUUxx * gamma_uu.xx;
	real gammaUUxxSq = gamma_uu.xx * gamma_uu.xx;
	
	real sqrt_f = eig.alpha_sqrt_f / eig.alpha;
	real f = sqrt_f * sqrt_f;
	// from 2004 Bona et al, "A symmetry breaking..." eqn A.20
	// mind you the 'm' in that form is for alpha_,t = -alpha^2 (f K - m Theta)
	// while in more modern Z4 papers it is alpha_,t = -alpha^2 f (K - m Theta)
	// the difference means that this eigendecomposition is probably incorrect.
	real lambda_1 = (2. - solver->m) / (f - 1);
	real lambda_2 = (2. * f - solver->m) / (f - 1);

	results.ptr[0] = (((lambda_2 * gamma_uu.xx * Z_l.x) - (lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.x.yy)) + ((lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.y.xy) - (2. * lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz)) + (lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.y.xz) + ((lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.z.xy) - (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.x.zz)) + (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.z.xz) + ((lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy)) + (lambda_2 * gamma_uu.xy * Z_l.y) + ((((2. * lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz)) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy)) - (lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz)) + ((lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz)) + (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz) + ((lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz)) + (lambda_2 * gamma_uu.xz * Z_l.z) + ((lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz) - (lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy)) + ((lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz)) + (sqrt_f * lambda_1 * Theta) + (sqrt_f * gamma_uu.xx * K_ll.xx) + (2. * sqrt_f * gamma_uu.xy * K_ll.xy) + (2. * sqrt_f * gamma_uu.xz * K_ll.xz) + (sqrt_f * gamma_uu.yy * K_ll.yy) + (2. * sqrt_f * gamma_uu.yz * K_ll.yz) + ((((sqrt_f * gamma_uu.zz * K_ll.zz) - (gamma_uu.xx * a_l.x)) - (gamma_uu.xy * a_l.y)) - (gamma_uu.xz * a_l.z)));
	results.ptr[1] = (Theta + ((gamma_uu.xx * Z_l.x) - (gamma_uu.xx * gamma_uu.yy * d_lll.x.yy)) + ((gamma_uu.xx * gamma_uu.yy * d_lll.y.xy) - (2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz)) + (gamma_uu.xx * gamma_uu.yz * d_lll.y.xz) + ((gamma_uu.xx * gamma_uu.yz * d_lll.z.xy) - (gamma_uu.xx * gamma_uu.zz * d_lll.x.zz)) + (gamma_uu.xx * gamma_uu.zz * d_lll.z.xz) + ((gamma_uu.xy * gamma_uu.xy * d_lll.x.yy) - (gamma_uu.xy * gamma_uu.xy * d_lll.y.xy)) + (gamma_uu.xy * Z_l.y) + ((((2. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz) - (gamma_uu.xy * gamma_uu.xz * d_lll.y.xz)) - (gamma_uu.xy * gamma_uu.xz * d_lll.z.xy)) - (gamma_uu.xy * gamma_uu.yz * d_lll.y.yz)) + ((gamma_uu.xy * gamma_uu.yz * d_lll.z.yy) - (gamma_uu.xy * gamma_uu.zz * d_lll.y.zz)) + (gamma_uu.xy * gamma_uu.zz * d_lll.z.yz) + ((gamma_uu.xz * gamma_uu.xz * d_lll.x.zz) - (gamma_uu.xz * gamma_uu.xz * d_lll.z.xz)) + (gamma_uu.xz * Z_l.z) + ((gamma_uu.xz * gamma_uu.yy * d_lll.y.yz) - (gamma_uu.xz * gamma_uu.yy * d_lll.z.yy)) + ((gamma_uu.xz * gamma_uu.yz * d_lll.y.zz) - (gamma_uu.xz * gamma_uu.yz * d_lll.z.yz)));
	results.ptr[2] = ((-(((a_l.y - (2. * K_ll.xy)) - (2. * Z_l.y)) + ((gamma_uu.xx * d_lll.y.xx) - (2. * gamma_uu.xy * d_lll.x.yy)) + ((2. * gamma_uu.xy * d_lll.y.xy) - (2. * gamma_uu.xz * d_lll.x.yz)) + (((2. * gamma_uu.xz * d_lll.y.xz) - (gamma_uu.yy * d_lll.y.yy)) - (2. * gamma_uu.yz * d_lll.z.yy)) + ((gamma_uu.zz * d_lll.y.zz) - (2. * gamma_uu.zz * d_lll.z.yz)))) / 2.);
	results.ptr[3] = ((-(((a_l.z - (2. * K_ll.xz)) - (2. * Z_l.z)) + ((gamma_uu.xx * d_lll.z.xx) - (2. * gamma_uu.xy * d_lll.x.yz)) + ((2. * gamma_uu.xy * d_lll.z.xy) - (2. * gamma_uu.xz * d_lll.x.zz)) + ((2. * gamma_uu.xz * d_lll.z.xz) - (2. * gamma_uu.yy * d_lll.y.yz)) + (((gamma_uu.yy * d_lll.z.yy) - (2. * gamma_uu.yz * d_lll.y.zz)) - (gamma_uu.zz * d_lll.z.zz)))) / 2.);
	results.ptr[4] = (((K_ll.yy - (gamma_uu.xx * d_lll.x.yy)) - (gamma_uu.xy * d_lll.y.yy)) - (gamma_uu.xz * d_lll.z.yy));
	results.ptr[5] = (((K_ll.yz - (gamma_uu.xx * d_lll.x.yz)) - (gamma_uu.xy * d_lll.y.yz)) - (gamma_uu.xz * d_lll.z.yz));
	results.ptr[6] = (((K_ll.zz - (gamma_uu.xx * d_lll.x.zz)) - (gamma_uu.xy * d_lll.y.zz)) - (gamma_uu.xz * d_lll.z.zz));
	results.ptr[7] = a_l.y;
	results.ptr[8] = a_l.z;
	results.ptr[9] = d_lll.y.xx;
	results.ptr[10] = d_lll.y.xy;
	results.ptr[11] = d_lll.y.xz;
	results.ptr[12] = d_lll.y.yy;
	results.ptr[13] = d_lll.y.yz;
	results.ptr[14] = d_lll.y.zz;
	results.ptr[15] = d_lll.z.xx;
	results.ptr[16] = d_lll.z.xy;
	results.ptr[17] = d_lll.z.xz;
	results.ptr[18] = d_lll.z.yy;
	results.ptr[19] = d_lll.z.yz;
	results.ptr[20] = d_lll.z.zz;
	results.ptr[21] = (a_l.x + ((Z_l.x - (f * gamma_uu.xx * d_lll.x.xx)) - (gamma_uu.xy * d_lll.x.xy)) + (((gamma_uu.xy * d_lll.y.xx) - (2. * gamma_uu.xy * f * d_lll.x.xy)) - (gamma_uu.xz * d_lll.x.xz)) + (((gamma_uu.xz * d_lll.z.xx) - (2. * gamma_uu.xz * f * d_lll.x.xz)) - (gamma_uu.yy * d_lll.x.yy)) + (((gamma_uu.yy * d_lll.y.xy) - (gamma_uu.yy * f * d_lll.x.yy)) - (2. * gamma_uu.yz * d_lll.x.yz)) + (gamma_uu.yz * d_lll.y.xz) + (((gamma_uu.yz * d_lll.z.xy) - (2. * gamma_uu.yz * f * d_lll.x.yz)) - (gamma_uu.zz * d_lll.x.zz)) + ((gamma_uu.zz * d_lll.z.xz) - (gamma_uu.zz * f * d_lll.x.zz)));
	results.ptr[22] = (a_l.y + (Z_l.y - (f * gamma_uu.yy * d_lll.y.yy)) + (((gamma_uu.xx * d_lll.x.xy) - (gamma_uu.xx * d_lll.y.xx)) - (gamma_uu.xx * f * d_lll.y.xx)) + (((gamma_uu.xy * d_lll.x.yy) - (gamma_uu.xy * d_lll.y.xy)) - (2. * gamma_uu.xy * f * d_lll.y.xy)) + ((gamma_uu.xz * d_lll.x.yz) - (2. * gamma_uu.xz * d_lll.y.xz)) + (((gamma_uu.xz * d_lll.z.xy) - (2. * gamma_uu.xz * f * d_lll.y.xz)) - (gamma_uu.yz * d_lll.y.yz)) + (((gamma_uu.yz * d_lll.z.yy) - (2. * gamma_uu.yz * f * d_lll.y.yz)) - (gamma_uu.zz * d_lll.y.zz)) + ((gamma_uu.zz * d_lll.z.yz) - (gamma_uu.zz * f * d_lll.y.zz)));
	results.ptr[23] = (a_l.z + (Z_l.z - (f * gamma_uu.zz * d_lll.z.zz)) + (((gamma_uu.xx * d_lll.x.xz) - (gamma_uu.xx * d_lll.z.xx)) - (gamma_uu.xx * f * d_lll.z.xx)) + (gamma_uu.xy * d_lll.x.yz) + (((gamma_uu.xy * d_lll.y.xz) - (2. * gamma_uu.xy * d_lll.z.xy)) - (2. * gamma_uu.xy * f * d_lll.z.xy)) + (((gamma_uu.xz * d_lll.x.zz) - (gamma_uu.xz * d_lll.z.xz)) - (2. * gamma_uu.xz * f * d_lll.z.xz)) + (((gamma_uu.yy * d_lll.y.yz) - (gamma_uu.yy * d_lll.z.yy)) - (gamma_uu.yy * f * d_lll.z.yy)) + (((gamma_uu.yz * d_lll.y.zz) - (gamma_uu.yz * d_lll.z.yz)) - (2. * gamma_uu.yz * f * d_lll.z.yz)));
	results.ptr[24] = ((a_l.y + ((2. * K_ll.xy) - (2. * Z_l.y)) + ((gamma_uu.xx * d_lll.y.xx) - (2. * gamma_uu.xy * d_lll.x.yy)) + ((2. * gamma_uu.xy * d_lll.y.xy) - (2. * gamma_uu.xz * d_lll.x.yz)) + (((2. * gamma_uu.xz * d_lll.y.xz) - (gamma_uu.yy * d_lll.y.yy)) - (2. * gamma_uu.yz * d_lll.z.yy)) + ((gamma_uu.zz * d_lll.y.zz) - (2. * gamma_uu.zz * d_lll.z.yz))) / 2.);
	results.ptr[25] = ((a_l.z + ((2. * K_ll.xz) - (2. * Z_l.z)) + ((gamma_uu.xx * d_lll.z.xx) - (2. * gamma_uu.xy * d_lll.x.yz)) + ((2. * gamma_uu.xy * d_lll.z.xy) - (2. * gamma_uu.xz * d_lll.x.zz)) + ((2. * gamma_uu.xz * d_lll.z.xz) - (2. * gamma_uu.yy * d_lll.y.yz)) + (((gamma_uu.yy * d_lll.z.yy) - (2. * gamma_uu.yz * d_lll.y.zz)) - (gamma_uu.zz * d_lll.z.zz))) / 2.);
	results.ptr[26] = (K_ll.yy + (gamma_uu.xx * d_lll.x.yy) + (gamma_uu.xy * d_lll.y.yy) + (gamma_uu.xz * d_lll.z.yy));
	results.ptr[27] = (K_ll.yz + (gamma_uu.xx * d_lll.x.yz) + (gamma_uu.xy * d_lll.y.yz) + (gamma_uu.xz * d_lll.z.yz));
	results.ptr[28] = (K_ll.zz + (gamma_uu.xx * d_lll.x.zz) + (gamma_uu.xy * d_lll.y.zz) + (gamma_uu.xz * d_lll.z.zz));
	results.ptr[29] = ((Theta - (gamma_uu.xx * Z_l.x)) + ((gamma_uu.xx * gamma_uu.yy * d_lll.x.yy) - (gamma_uu.xx * gamma_uu.yy * d_lll.y.xy)) + (((2. * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz) - (gamma_uu.xx * gamma_uu.yz * d_lll.y.xz)) - (gamma_uu.xx * gamma_uu.yz * d_lll.z.xy)) + (((gamma_uu.xx * gamma_uu.zz * d_lll.x.zz) - (gamma_uu.xx * gamma_uu.zz * d_lll.z.xz)) - (gamma_uu.xy * gamma_uu.xy * d_lll.x.yy)) + (((gamma_uu.xy * gamma_uu.xy * d_lll.y.xy) - (gamma_uu.xy * Z_l.y)) - (2. * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz)) + (gamma_uu.xy * gamma_uu.xz * d_lll.y.xz) + (gamma_uu.xy * gamma_uu.xz * d_lll.z.xy) + ((gamma_uu.xy * gamma_uu.yz * d_lll.y.yz) - (gamma_uu.xy * gamma_uu.yz * d_lll.z.yy)) + (((gamma_uu.xy * gamma_uu.zz * d_lll.y.zz) - (gamma_uu.xy * gamma_uu.zz * d_lll.z.yz)) - (gamma_uu.xz * gamma_uu.xz * d_lll.x.zz)) + (((gamma_uu.xz * gamma_uu.xz * d_lll.z.xz) - (gamma_uu.xz * Z_l.z)) - (gamma_uu.xz * gamma_uu.yy * d_lll.y.yz)) + ((gamma_uu.xz * gamma_uu.yy * d_lll.z.yy) - (gamma_uu.xz * gamma_uu.yz * d_lll.y.zz)) + (gamma_uu.xz * gamma_uu.yz * d_lll.z.yz));
	results.ptr[30] = (-(((lambda_2 * gamma_uu.xx * Z_l.x) - (lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.x.yy)) + ((lambda_2 * gamma_uu.xx * gamma_uu.yy * d_lll.y.xy) - (2. * lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.x.yz)) + (lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.y.xz) + ((lambda_2 * gamma_uu.xx * gamma_uu.yz * d_lll.z.xy) - (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.x.zz)) + (lambda_2 * gamma_uu.xx * gamma_uu.zz * d_lll.z.xz) + ((lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.x.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.xy * d_lll.y.xy)) + (lambda_2 * gamma_uu.xy * Z_l.y) + ((((2. * lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.x.yz) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.y.xz)) - (lambda_2 * gamma_uu.xy * gamma_uu.xz * d_lll.z.xy)) - (lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.y.yz)) + ((lambda_2 * gamma_uu.xy * gamma_uu.yz * d_lll.z.yy) - (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.y.zz)) + (lambda_2 * gamma_uu.xy * gamma_uu.zz * d_lll.z.yz) + ((lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.x.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.xz * d_lll.z.xz)) + (lambda_2 * gamma_uu.xz * Z_l.z) + ((lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.y.yz) - (lambda_2 * gamma_uu.xz * gamma_uu.yy * d_lll.z.yy)) + ((((((((((((lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.y.zz) - (lambda_2 * gamma_uu.xz * gamma_uu.yz * d_lll.z.yz)) - (sqrt_f * lambda_1 * Theta)) - (sqrt_f * gamma_uu.xx * K_ll.xx)) - (2. * sqrt_f * gamma_uu.xy * K_ll.xy)) - (2. * sqrt_f * gamma_uu.xz * K_ll.xz)) - (sqrt_f * gamma_uu.yy * K_ll.yy)) - (2. * sqrt_f * gamma_uu.yz * K_ll.yz)) - (sqrt_f * gamma_uu.zz * K_ll.zz)) - (gamma_uu.xx * a_l.x)) - (gamma_uu.xy * a_l.y)) - (gamma_uu.xz * a_l.z))));

#endif
	return results;
}

//TODO these were based no noZeroRowsInFlux==false (I think) so maybe/certainly they are out of date
<?=eqn.cons_t?> eigen_rightTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.waves_t?> input,
	real3 x,
	normal_t n
) {
	<?=eqn.cons_t?> resultU = (<?=eqn.cons_t?>){.ptr={0. / 0.}};
	//for (int j = 0; j < numStates; ++j) {
	//	resultU.ptr[j] = 0;
	//}

#error here
#if 0	
	sym3 gamma_ll = sym3_swap(eig.gamma_ll, n.side);
	sym3 gamma_uu = sym3_swap(eig.gamma_uu, n.side);
	
	real sqrt_gammaUUxx = sqrt(gamma_uu.xx);
	real gammaUUxx_toThe_3_2 = sqrt_gammaUUxx * gamma_uu.xx;
	real gammaUUxxSq = gamma_uu.xx * gamma_uu.xx;

	real sqrt_f = eig.alpha_sqrt_f / eig.alpha;
	real f = sqrt_f * sqrt_f;
	real fSq = f * f;
	// from 2004 Bona et al, "A symmetry breaking..." eqn A.20
	real lambda_1 = (2. - solver->m) / (f - 1);
	real lambda_2 = (2. * f - solver->m) / (f - 1);
	
	resultU.ptr[0] = ((((input.ptr[0] - input.ptr[30]) - (lambda_2 * input.ptr[1])) + (lambda_2 * input.ptr[29]) + (2. * gamma_uu.xy * input.ptr[7]) + (2. * gamma_uu.xz * input.ptr[8])) / (-(2. * gamma_uu.xx)));
	resultU.ptr[1] = input.ptr[7];
	resultU.ptr[2] = input.ptr[8];
	resultU.ptr[3] = ((((2. * input.ptr[0]) - (2. * input.ptr[1])) + (4. * input.ptr[21] * gamma_uu.xx) + (((2. * input.ptr[29]) - (2. * input.ptr[30])) - (2. * lambda_2 * input.ptr[1])) + (((2. * lambda_2 * input.ptr[29]) - (4. * gamma_uu.xy * input.ptr[2] * f)) - (12. * gamma_uu.xy * input.ptr[7] * f)) + (8. * gamma_uu.xy * input.ptr[9] * gamma_uu.xx * f) + (8. * gamma_uu.xy * gamma_uu.xy * input.ptr[10] * f) + (4. * gamma_uu.xy * input.ptr[22]) + (4. * gamma_uu.xy * input.ptr[24] * f) + (8. * gamma_uu.xy * fSq * input.ptr[9] * gamma_uu.xx) + (16. * gamma_uu.xy * gamma_uu.xy * fSq * input.ptr[10]) + (8. * gamma_uu.xy * f * input.ptr[22]) + (8. * gamma_uu.xy * gamma_uu.xz * input.ptr[11] * f) + (8. * gamma_uu.xy * gamma_uu.xz * input.ptr[16] * f) + (16. * gamma_uu.xy * gamma_uu.xz * fSq * input.ptr[11]) + (16. * gamma_uu.xy * gamma_uu.xz * fSq * input.ptr[16]) + (8. * gamma_uu.xy * gamma_uu.yz * input.ptr[13] * f) + (((16. * gamma_uu.xy * gamma_uu.yz * fSq * input.ptr[13]) - (4. * gamma_uu.xz * input.ptr[3] * f)) - (12. * gamma_uu.xz * input.ptr[8] * f)) + (8. * gamma_uu.xz * input.ptr[15] * gamma_uu.xx * f) + (8. * gamma_uu.xz * gamma_uu.xz * input.ptr[17] * f) + (4. * gamma_uu.xz * input.ptr[23]) + (4. * gamma_uu.xz * input.ptr[25] * f) + (8. * gamma_uu.xz * fSq * input.ptr[15] * gamma_uu.xx) + (16. * gamma_uu.xz * gamma_uu.xz * fSq * input.ptr[17]) + (8. * gamma_uu.xz * f * input.ptr[23]) + (8. * gamma_uu.xz * gamma_uu.yz * input.ptr[19] * f) + ((16. * gamma_uu.xz * gamma_uu.yz * fSq * input.ptr[19]) - (2. * gamma_uu.yy * input.ptr[4] * f)) + (2. * gamma_uu.yy * input.ptr[26] * f) + (4. * gamma_uu.yy * gamma_uu.xy * input.ptr[12] * f) + (8. * gamma_uu.yy * gamma_uu.xy * fSq * input.ptr[12]) + (4. * gamma_uu.yy * gamma_uu.xz * input.ptr[18] * f) + ((8. * gamma_uu.yy * gamma_uu.xz * fSq * input.ptr[18]) - (4. * gamma_uu.yz * input.ptr[5] * f)) + ((4. * gamma_uu.yz * input.ptr[27] * f) - (2. * gamma_uu.zz * input.ptr[6] * f)) + (2. * gamma_uu.zz * input.ptr[28] * f) + (4. * gamma_uu.zz * gamma_uu.xy * input.ptr[14] * f) + (8. * gamma_uu.zz * gamma_uu.xy * fSq * input.ptr[14]) + (4. * gamma_uu.zz * gamma_uu.xz * input.ptr[20] * f) + (8. * gamma_uu.zz * gamma_uu.xz * fSq * input.ptr[20])) / (-(4. * gammaUUxxSq * f)));
	resultU.ptr[4] = ((-(input.ptr[2] + (((((((3. * input.ptr[7]) - (input.ptr[9] * gamma_uu.xx)) - (2. * input.ptr[22])) - input.ptr[24]) - (2. * f * input.ptr[9] * gamma_uu.xx)) - (4. * gamma_uu.xy * f * input.ptr[10])) - (2. * gamma_uu.xz * input.ptr[11])) + ((((((((2. * gamma_uu.xz * input.ptr[16]) - (4. * gamma_uu.xz * f * input.ptr[11])) - (gamma_uu.yy * input.ptr[12])) - (2. * gamma_uu.yy * f * input.ptr[12])) - (2. * gamma_uu.yz * input.ptr[13])) - (4. * gamma_uu.yz * f * input.ptr[13])) - (gamma_uu.zz * input.ptr[14])) - (2. * gamma_uu.zz * f * input.ptr[14])))) / (2. * gamma_uu.xx));
	resultU.ptr[5] = ((-(input.ptr[3] + (((((3. * input.ptr[8]) - (input.ptr[15] * gamma_uu.xx)) - (2. * input.ptr[23])) - input.ptr[25]) - (2. * f * input.ptr[15] * gamma_uu.xx)) + ((((((((((2. * gamma_uu.xy * input.ptr[11]) - (2. * gamma_uu.xy * input.ptr[16])) - (4. * gamma_uu.xy * f * input.ptr[16])) - (4. * gamma_uu.xz * f * input.ptr[17])) - (gamma_uu.yy * input.ptr[18])) - (2. * gamma_uu.yy * f * input.ptr[18])) - (2. * gamma_uu.yz * input.ptr[19])) - (4. * gamma_uu.yz * f * input.ptr[19])) - (gamma_uu.zz * input.ptr[20])) - (2. * gamma_uu.zz * f * input.ptr[20])))) / (2. * gamma_uu.xx));
	resultU.ptr[6] = ((-((input.ptr[4] - input.ptr[26]) + (2. * gamma_uu.xy * input.ptr[12]) + (2. * gamma_uu.xz * input.ptr[18]))) / (2. * gamma_uu.xx));
	resultU.ptr[7] = ((-((input.ptr[5] - input.ptr[27]) + (2. * gamma_uu.xy * input.ptr[13]) + (2. * gamma_uu.xz * input.ptr[19]))) / (2. * gamma_uu.xx));
	resultU.ptr[8] = ((-((input.ptr[6] - input.ptr[28]) + (2. * gamma_uu.xy * input.ptr[14]) + (2. * gamma_uu.xz * input.ptr[20]))) / (2. * gamma_uu.xx));
	resultU.ptr[9] = input.ptr[9];
	resultU.ptr[10] = input.ptr[10];
	resultU.ptr[11] = input.ptr[11];
	resultU.ptr[12] = input.ptr[12];
	resultU.ptr[13] = input.ptr[13];
	resultU.ptr[14] = input.ptr[14];
	resultU.ptr[15] = input.ptr[15];
	resultU.ptr[16] = input.ptr[16];
	resultU.ptr[17] = input.ptr[17];
	resultU.ptr[18] = input.ptr[18];
	resultU.ptr[19] = input.ptr[19];
	resultU.ptr[20] = input.ptr[20];
	resultU.ptr[21] = ((input.ptr[0] + ((((((((((((input.ptr[30] - (lambda_1 * input.ptr[1] * sqrt_f)) - (lambda_1 * input.ptr[29] * sqrt_f)) - (2. * gamma_uu.xy * input.ptr[2] * sqrt_f)) - (2. * gamma_uu.xy * input.ptr[24] * sqrt_f)) - (2. * gamma_uu.xz * input.ptr[3] * sqrt_f)) - (2. * gamma_uu.xz * input.ptr[25] * sqrt_f)) - (gamma_uu.yy * input.ptr[4] * sqrt_f)) - (gamma_uu.yy * input.ptr[26] * sqrt_f)) - (2. * gamma_uu.yz * input.ptr[5] * sqrt_f)) - (2. * gamma_uu.yz * input.ptr[27] * sqrt_f)) - (gamma_uu.zz * input.ptr[6] * sqrt_f)) - (gamma_uu.zz * input.ptr[28] * sqrt_f))) / (2. * sqrt_f * gamma_uu.xx));
	resultU.ptr[22] = ((input.ptr[2] + input.ptr[24]) / 2.);
	resultU.ptr[23] = ((input.ptr[3] + input.ptr[25]) / 2.);
	resultU.ptr[24] = ((input.ptr[4] + input.ptr[26]) / 2.);
	resultU.ptr[25] = ((input.ptr[5] + input.ptr[27]) / 2.);
	resultU.ptr[26] = ((input.ptr[6] + input.ptr[28]) / 2.);
	resultU.ptr[27] = ((input.ptr[1] + input.ptr[29]) / 2.);
	resultU.ptr[28] = ((((((input.ptr[1] - input.ptr[29]) - (gamma_uu.xy * input.ptr[2])) - (gamma_uu.xy * input.ptr[7])) - (gamma_uu.xy * input.ptr[9] * gamma_uu.xx)) + (((((gamma_uu.xy * input.ptr[24]) - (2. * gamma_uu.xy * gamma_uu.yz * input.ptr[13])) - (gamma_uu.xz * input.ptr[3])) - (gamma_uu.xz * input.ptr[8])) - (gamma_uu.xz * input.ptr[15] * gamma_uu.xx)) + (((gamma_uu.xz * input.ptr[25]) - (gamma_uu.yy * input.ptr[4])) - (2. * gamma_uu.yy * input.ptr[10] * gamma_uu.xx)) + ((((((gamma_uu.yy * input.ptr[26]) - (gamma_uu.yy * gamma_uu.xy * input.ptr[12])) - (gamma_uu.yy * gamma_uu.xz * input.ptr[18])) - (2. * gamma_uu.yz * input.ptr[5])) - (2. * gamma_uu.yz * input.ptr[11] * gamma_uu.xx)) - (2. * gamma_uu.yz * input.ptr[16] * gamma_uu.xx)) + ((((2. * gamma_uu.yz * input.ptr[27]) - (2. * gamma_uu.yz * gamma_uu.xz * input.ptr[19])) - (gamma_uu.zz * input.ptr[6])) - (2. * gamma_uu.zz * input.ptr[17] * gamma_uu.xx)) + (((gamma_uu.zz * input.ptr[28]) - (gamma_uu.zz * gamma_uu.xy * input.ptr[14])) - (gamma_uu.zz * gamma_uu.xz * input.ptr[20]))) / (2. * gamma_uu.xx));
	resultU.ptr[29] = (((input.ptr[2] * gamma_uu.xx) + ((input.ptr[7] * gamma_uu.xx) - (input.ptr[24] * gamma_uu.xx)) + (((gammaUUxxSq * input.ptr[9]) - (gamma_uu.xx * gamma_uu.yy * input.ptr[12])) - (2. * gamma_uu.xx * gamma_uu.yz * input.ptr[18])) + (gamma_uu.xy * input.ptr[4]) + (2. * gamma_uu.xy * input.ptr[10] * gamma_uu.xx) + ((2. * gamma_uu.xy * gamma_uu.xy * input.ptr[12]) - (gamma_uu.xy * input.ptr[26])) + (2. * gamma_uu.xy * gamma_uu.xz * input.ptr[13]) + (gamma_uu.xz * input.ptr[5]) + (2. * gamma_uu.xz * input.ptr[11] * gamma_uu.xx) + ((2. * gamma_uu.xz * gamma_uu.xz * input.ptr[19]) - (gamma_uu.xz * input.ptr[27])) + (2. * gamma_uu.xz * gamma_uu.xy * input.ptr[18]) + ((gamma_uu.zz * input.ptr[14] * gamma_uu.xx) - (2. * gamma_uu.zz * gamma_uu.xx * input.ptr[19]))) / (2. * gamma_uu.xx));
	resultU.ptr[30] = (((input.ptr[3] * gamma_uu.xx) + ((input.ptr[8] * gamma_uu.xx) - (input.ptr[25] * gamma_uu.xx)) + ((((gammaUUxxSq * input.ptr[15]) - (2. * gamma_uu.xx * gamma_uu.yy * input.ptr[13])) - (2. * gamma_uu.xx * gamma_uu.yz * input.ptr[14])) - (gamma_uu.xx * gamma_uu.zz * input.ptr[20])) + (gamma_uu.xy * input.ptr[5]) + (2. * gamma_uu.xy * gamma_uu.xy * input.ptr[13]) + ((2. * gamma_uu.xy * input.ptr[16] * gamma_uu.xx) - (gamma_uu.xy * input.ptr[27])) + (2. * gamma_uu.xy * gamma_uu.xz * input.ptr[19]) + (gamma_uu.xz * input.ptr[6]) + (2. * gamma_uu.xz * input.ptr[17] * gamma_uu.xx) + ((2. * gamma_uu.xz * gamma_uu.xz * input.ptr[20]) - (gamma_uu.xz * input.ptr[28])) + (2. * gamma_uu.xz * gamma_uu.xy * input.ptr[14]) + (gamma_uu.yy * input.ptr[18] * gamma_uu.xx)) / (2. * gamma_uu.xx));

	resultU.a_l = real3_swap(resultU.a_l, n.side);			//0-2
	resultU.d_lll = _3sym3_swap(resultU.d_lll, n.side);		//3-20
	resultU.K_ll = sym3_swap(resultU.K_ll, n.side);			//21-26
	resultU.Theta = resultU.Theta;							//27
	resultU.Z_l = real3_swap(resultU.Z_l, n.side);			//28-30
#endif
	return resultU;
}
]],
	}

	-- used by roe, some plm
	solver.modules:add{
		name = 'eigen_fluxTransform',
		depends = {
			'solver.solver_t',
			'eqn.cons_t',
			'normal_t',
		},
		code = self:template[[
//so long as roeUseFluxFromCons isn't set for the roe solver, 
// and fluxFromCons is provided/unused,
// eigen_fluxTransform isn't needed.
// but some solvers do use a boilerplate right(lambda(left(U)))
//however if you want to use the HLL solver then fluxFromCons is needed
//...however fluxFromCons is not provided by this eqn.
<?=eqn.cons_t?> eigen_fluxTransform(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.eigen_t?> eig,
	<?=eqn.cons_t?> inputU,
	real3 x,
	normal_t n
) {
#if 0	
	//default
	<?=eqn.waves_t?> waves = eigen_leftTransform(solver, eig, inputU, x, n);
	<?=eqn:eigenWaveCodePrefix('n', 'eig', 'x')?>
<? for j=0,eqn.numWaves-1 do 
?>	waves.ptr[<?=j?>] *= <?=eqn:eigenWaveCode('n', 'eig', 'x', j)?>;
<? end 
?>	return eigen_rightTransform(solver, eig, waves, x, n);
#else
	<?=eqn.cons_t?> F = {.ptr={0. / 0.}};
	return F;
#endif
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
		real sqrt_gammaUjj = sqrt(gammaUjj);
		real lambdaLight = sqrt_gammaUjj * U->alpha;
		real lambdaGauge = sqrt_gammaUjj * alpha_sqrt_f;
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

<? if solver.coord.vectorComponent == 'cartesian' then ?>

//taken from sym3 parallel propagate
//so TODO somehow generalize this index transform for all data types
//but this looks like a job for lambda functions, not so much C ... maybe I can generate it with scripts ... 

sym3 sym3_rotateFrom(sym3 m, real3 n) {
	real3x3 t = real3x3_from_sym3(m);
	t.x = real3_rotateFrom(t.x, n);
	t.y = real3_rotateFrom(t.y, n);
	t.z = real3_rotateFrom(t.z, n);
	t = real3x3_transpose(t);
	t.x = real3_rotateFrom(t.x, n);
	t.y = real3_rotateFrom(t.y, n);
	t.z = real3_rotateFrom(t.z, n);
	return sym3_from_real3x3(t);
}

sym3 sym3_rotateTo(sym3 m, real3 n) {
	real3x3 t = real3x3_from_sym3(m);
	t.x = real3_rotateTo(t.x, n);
	t.y = real3_rotateTo(t.y, n);
	t.z = real3_rotateTo(t.z, n);
	t = real3x3_transpose(t);
	t.x = real3_rotateTo(t.x, n);
	t.y = real3_rotateTo(t.y, n);
	t.z = real3_rotateTo(t.z, n);
	return sym3_from_real3x3(t);
}

_3sym3 _3sym3_rotateFrom(_3sym3 m, real3 n) {
	real3x3x3 t = real3x3x3_from__3sym3(m);
	real3 tmp;
<?	local is = require 'ext.table'()
	for e=1,3 do
		for i,xi in ipairs(xNames) do
			is[e] = xi
			for j,xj in ipairs(xNames) do
				is[e%3+1] = xj
				for k,xk in ipairs(xNames) do
					is[(e+1)%3+1] = xk
?>	tmp.<?=xk?> = t.<?=is:concat'.'?>;
<?				end
?>	tmp = real3_rotateFrom(tmp, n);
<?				for k,xk in ipairs(xNames) do
					is[(e+1)%3+1] = xk
?>	t.<?=is:concat'.'?> = tmp.<?=xk?>;
<?				end
			end
		end
	end
?>	return _3sym3_from_real3x3x3(t);
}

_3sym3 _3sym3_rotateTo(_3sym3 m, real3 n) {
	real3x3x3 t = real3x3x3_from__3sym3(m);
	real3 tmp;
<?	local is = require 'ext.table'()
	for e=1,3 do
		for i,xi in ipairs(xNames) do
			is[e] = xi
			for j,xj in ipairs(xNames) do
				is[e%3+1] = xj
				for k,xk in ipairs(xNames) do
					is[(e+1)%3+1] = xk
?>	tmp.<?=xk?> = t.<?=is:concat'.'?>;
<?				end
?>	tmp = real3_rotateTo(tmp, n);
<?				for k,xk in ipairs(xNames) do
					is[(e+1)%3+1] = xk
?>	t.<?=is:concat'.'?> = tmp.<?=xk?>;
<?				end
			end
		end
	end
?>	return _3sym3_from_real3x3x3(t);
}

<? end ?>

<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	real f_alpha = calc_f_alpha(U.alpha);
	
	real det_gamma = sym3_det(U.gamma_ll);
	
	real alpha = U.alpha;
	real Theta = U.Theta;
	
	real3 Z_l = U.Z_l;
	real3 a_l = U.a_l;
	sym3 gamma_ll = U.gamma_ll;
	sym3 K_ll = U.K_ll;
	_3sym3 d_lll = U.d_lll;

<? if solver.coord.vectorComponent == 'cartesian' then ?>

// I'm using .side for holonomic(coordinate) and anholonomic(orthonormal)
//but for cartesian vector componets there is no .side, just .n, which is covariant iirc
//and I haven't derived the flux in arbitrary-normal form, just in x-axis form (and i swap x<->y or z to calculate their fluxes)
//so here I'm going to just wing it
	real3 n_l = normal_l1(n);

	Z_l = real3_rotateFrom(Z_l, n_l);
	a_l = real3_rotateFrom(a_l, n_l);
	gamma_ll = sym3_rotateFrom(gamma_ll, n_l);
	K_ll = sym3_rotateFrom(K_ll, n_l);
	d_lll = _3sym3_rotateFrom(d_lll, n_l);

<? else ?>

	if (false) {}
	<? for side=0,solver.dim-1 do ?>
	else if (n.side == <?=side?>) {
		Z_l = real3_swap<?=side?>(Z_l);
		a_l = real3_swap<?=side?>(a_l);
		gamma_ll = sym3_swap<?=side?>(gamma_ll);
		K_ll = sym3_swap<?=side?>(K_ll);
		d_lll = _3sym3_swap<?=side?>(d_lll);
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
	
<? end ?>

	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma);

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

<? if solver.coord.vectorComponent == 'cartesian' then ?>

	F.Z_l = real3_rotateFrom(F.Z_l, n_l);
	F.a_l = real3_rotateFrom(F.a_l, n_l);
	F.gamma_ll = sym3_rotateFrom(F.gamma_ll, n_l);
	F.K_ll = sym3_rotateFrom(F.K_ll, n_l);
	F.d_lll = _3sym3_rotateFrom(F.d_lll, n_l);
	
<? else ?>

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

<? end ?>

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
		'sym3sym3',	-- d_kij,l
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
	*U = (<?=eqn.cons_t?>){.ptr={ 0. / 0. }};

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

	*U = (<?=eqn.cons_t?>){.ptr={ 0. / 0. }};
	
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

	vars:insert{name='alpha vs a_i', code=self:template[[
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
]], type='real3'}

	-- d_kij = gamma_ij,k
	for i,xi in ipairs(xNames) do
		vars:insert{name='gamma_ij vs d_'..xi..'ij', code=self:template([[
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
}), type='sym3'}
	end

	return vars
end

function Z4_2004Bona:eigenWaveCodePrefix(n, eig, x, waveIndex)
	return self:template([[
	real eig_lambdaLight = <?=eig?>.sqrt_gammaUnn * <?=eig?>.alpha;
	real eig_lambdaGauge = <?=eig?>.sqrt_gammaUnn * <?=eig?>.alpha_sqrt_f;
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
	return 'min(0., min(-'..betaUi..' - eig_lambdaGauge, -'..betaUi..' - eig_lambdaLight))'
end

function Z4_2004Bona:eigenWaveMaxCode(n, eig, x)
	local betaUi
	if self.useShift ~= 'none' then
		betaUi = '('..eig..').beta_u.s['..n..'.side]'
	else
		betaUi = '0'
	end
	return 'max(0., max(-'..betaUi..' + eig_lambdaGauge, -'..betaUi..' + eig_lambdaLight))'
end

function Z4_2004Bona:consWaveCodePrefix(n, U, x, waveIndex)
	return self:template([[
	real det_gamma = sym3_det(<?=U?>.gamma_ll);
	sym3 gamma_uu = sym3_inv(<?=U?>.gamma_ll, det_gamma);
	
<? if solver.coord.vectorComponent == 'cartesian' then ?>
	real3 n_l = normal_l1(n);
	real gammaUnn = real3_weightedLenSq(n_l, gamma_uu);
<? else ?>
	real gammaUnn = 0./0.;
	if (n.side == 0) {
		gammaUnn = gamma_uu.xx;
	} else if (n.side == 1) {
		gammaUnn = gamma_uu.yy;
	} else if (n.side == 2) {
		gammaUnn = gamma_uu.zz;
	}
<? end ?>

	real sqrt_gammaUnn = sqrt(gammaUnn);
	real eig_lambdaLight = sqrt_gammaUnn * <?=U?>.alpha;
	real alpha_sqrt_f = sqrt(calc_f_alphaSq(<?=U?>.alpha));
	real eig_lambdaGauge = sqrt_gammaUnn * alpha_sqrt_f;
]], {
		U = '('..U..')',
		n = n,
	})
end
Z4_2004Bona.consWaveCode = Z4_2004Bona.eigenWaveCode

return Z4_2004Bona

local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local makestruct = require 'eqn.makestruct'
local Equation = require 'eqn.eqn'

local GLM_MHD = class(Equation)

GLM_MHD.name = 'GLM-MHD'

GLM_MHD.primVars = table{
	{rho = 'real'},
	{v = 'real3'},
	{P = 'real'},
	{B = 'real3'},
	{psi = 'real'},
}

GLM_MHD.consVars = table{
	{rho = 'real'},
	{m = 'real3'},
	{ETotal = 'real'},
	{B = 'real3'},
	{psi = 'real'},
}

GLM_MHD.mirrorVars = {{'m.x', 'B.x'}, {'m.y', 'B.y'}, {'m.z', 'B.z'}}

GLM_MHD.hasEigenCode = true
GLM_MHD.hasFluxFromConsCode = true
GLM_MHD.roeUseFluxFromCons = true
GLM_MHD.useSourceTerm = true
GLM_MHD.useConstrianU = true

GLM_MHD.useFixedCh = false	-- true = use a gui var, false = calculate by max(|v_i|+Cf)

-- hmm, we want init.euler and init.mhd here ...
GLM_MHD.initStates = require 'init.euler'

function GLM_MHD:init(args)
	GLM_MHD.super.init(self, args)

	local UpdatePsi = require 'op.glm-mhd-update-psi'
	self.solver.ops:insert(UpdatePsi{solver=self.solver})
end

GLM_MHD.guiVars = table{
	{name='heatCapacityRatio', value=2},	-- 5/3 for most problems, but 2 for Brio-Wu, so I will just set it here for now (in case something else is reading it before it is set there)
	{name='mu0', value=1},	-- this should be 4 pi for natural units, but I haven't verified that all mu0's are where they should be ...
	{name='Cp', value=1, compileTime=true},
}

if GLM_MHD.useFixedCh then
	GLM_MHD.guiVars:insert{name='Ch', value=.1}
end

function GLM_MHD:getCommonFuncCode()
	return template([[
static inline real calc_eKin(<?=eqn.prim_t?> W) { return .5 * real3_lenSq(W.v); }
static inline real calc_EKin(<?=eqn.prim_t?> W) { return W.rho * calc_eKin(W); }
static inline real calc_EInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return W.P / (solver->heatCapacityRatio - 1.); }
static inline real calc_eInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_EInt(solver, W) / W.rho; }
static inline real calc_EMag(<?=eqn.prim_t?> W) { return .5 * real3_lenSq(W.B); }
static inline real calc_eMag(<?=eqn.prim_t?> W) { return calc_EMag(W) / W.rho; }
static inline real calc_PMag(<?=eqn.prim_t?> W) { return .5 * real3_lenSq(W.B); }
static inline real calc_EHydro(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_EKin(W) + calc_EInt(solver, W); }
static inline real calc_eHydro(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_EHydro(solver, W) / W.rho; }
static inline real calc_ETotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_EKin(W) + calc_EInt(solver, W) + calc_EMag(W); }
static inline real calc_eTotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_ETotal(solver, W) / W.rho; }
static inline real calc_H(constant <?=solver.solver_t?>* solver, real P) { return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)); }
static inline real calc_h(constant <?=solver.solver_t?>* solver, real rho, real P) { return calc_H(solver, P) / rho; }
static inline real calc_HTotal(<?=eqn.prim_t?> W, real ETotal) { return W.P + calc_PMag(W) + ETotal; }
static inline real calc_hTotal(<?=eqn.prim_t?> W, real ETotal) { return calc_HTotal(W, ETotal) / W.rho; }
static inline real calc_Cs(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return sqrt(solver->heatCapacityRatio * W.P / W.rho); }
]], {
		solver = self.solver,
		eqn = self,
	})
end

function GLM_MHD:getPrimConsCode()
	return template([[
static inline <?=eqn.prim_t?> primFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W;
	W.rho = U.rho;
	W.v = real3_real_mul(U.m, 1./U.rho);
	W.B = U.B;
	real vSq = real3_lenSq(W.v);
	real BSq = real3_lenSq(W.B);
	real EKin = .5 * U.rho * vSq;
	real EMag = .5 * BSq;
	real EInt = U.ETotal - EKin - EMag;
	W.P = EInt * (solver->heatCapacityRatio - 1.);
	W.P = max(W.P, (real)1e-7);
	W.rho = max(W.rho, (real)1e-7);
	W.psi = U.psi;
	return W;
}

static inline <?=eqn.cons_t?> consFromPrim(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> W,
	real3 x
) {
	<?=eqn.cons_t?> U;
	U.rho = W.rho;
	U.m = real3_real_mul(W.v, W.rho);
	U.B = W.B;
	real vSq = real3_lenSq(W.v);
	real BSq = real3_lenSq(W.B);
	real EKin = .5 * W.rho * vSq;
	real EMag = .5 * BSq;
	real EInt = W.P / (solver->heatCapacityRatio - 1.);
	U.ETotal = EInt + EKin + EMag;
	U.psi = W.psi;
	return U;
}

<?=eqn.cons_t?> apply_dU_dW(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA, 
	<?=eqn.prim_t?> W, 
	real3 x
) {
	return (<?=eqn.cons_t?>){
		.rho = W.rho,
		.m = real3_add(
			real3_real_mul(WA.v, W.rho),
			real3_real_mul(W.v, WA.rho)),
		.B = WA.B,
		.ETotal = W.rho * .5 * real3_dot(WA.v, WA.v)
			+ WA.rho * real3_dot(W.v, WA.v)
			+ real3_dot(W.B, WA.B) / solver->mu0
			+ W.P / (solver->heatCapacityRatio - 1.),
		.psi = W.psi,
	};
}

<?=eqn.prim_t?> apply_dW_dU(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA,
	<?=eqn.cons_t?> U,
	real3 x
) {
	return (<?=eqn.prim_t?>){
		.rho = U.rho,
		.v = real3_sub(
			real3_real_mul(U.m, 1. / WA.rho),
			real3_real_mul(WA.v, U.rho / WA.rho)),
		.B = U.B,
		.P = (solver->heatCapacityRatio - 1.) *  (
			.5 * U.rho * real3_dot(WA.v, WA.v)
			- real3_dot(U.m, WA.v)
			- real3_dot(U.B, WA.B) / solver->mu0
			+ U.ETotal),
		.psi = U.psi,
	};
}
]], {
	solver = self.solver,
	eqn = self,
})
end

GLM_MHD.initStateCode = [[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
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
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real3 B = real3_zero;

	<?=code?>
	
	<?=eqn.prim_t?> W = {.rho=rho, .v=v, .P=P, .B=B, .psi=0};
	UBuf[index] = consFromPrim(solver, W, x);
}

kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;

	U->psi = .5 * (0.
<? 
for j=0,solver.dim-1 do 
?>		+ (U[solver->stepsize.s<?=j?>].B.s<?=j?> 
			- U[-solver->stepsize.s<?=j?>].B.s<?=j?>
		) / solver->grid_dx.s<?=j?>
<? 
end 
?>	);
}
]]

GLM_MHD.solverCodeFile = 'eqn/glm-mhd.cl'

GLM_MHD.displayVarCodeUsesPrims = true

function GLM_MHD:getDisplayVars()
	return GLM_MHD.super.getDisplayVars(self):append{
		{v = '*value_real3 = W.v;', type='real3'},
		{['div B'] = template([[
	*value = .5 * (0.
<? 
for j=0,solver.dim-1 do 
?>		+ (U[solver->stepsize.s<?=j?>].<?=field?>.s<?=j?> 
			- U[-solver->stepsize.s<?=j?>].<?=field?>.s<?=j?>
		) / solver->grid_dx.s<?=j?>
<? 
end 
?>	);
]], {solver=self.solver, field='B'})},
		{P = '*value = W.P;'},
		--{PMag = '*value = calc_PMag(W);'},
		--{PTotal = '*value = W.P + calc_PMag(W);'},
		--{eInt = '*value = calc_eInt(solver, W);'},
		{EInt = '*value = calc_EInt(solver, W);'},
		--{eKin = '*value = calc_eKin(W);'},
		{EKin = '*value = calc_EKin(W);'},
		--{eHydro = '*value = calc_eHydro(solver, W);'},
		{EHydro = '*value = calc_EHydro(solver, W);'},
		--{eMag = '*value = calc_eMag(W);'},
		{EMag = '*value = calc_EMag(W);'},
		--{eTotal = '*value = U->ETotal / W.rho;'},
		{S = '*value = W.P / pow(W.rho, (real)solver->heatCapacityRatio);'},
		{H = '*value = calc_H(solver, W.P);'},
		--{h = '*value = calc_H(solver, W.P) / W.rho;'},
		--{HTotal = '*value = calc_HTotal(W, U->ETotal);'},
		--{hTotal = '*value = calc_hTotal(W, U->ETotal);'},
		--{Cs = '*value = calc_Cs(solver, W); },
		{['primitive reconstruction error'] = template([[
		//prim have just been reconstructed from cons
		//so reconstruct cons from prims again and calculate the difference
		<?=eqn.cons_t?> U2 = consFromPrim(solver, W, x);
		*value = 0;
		for (int j = 0; j < numIntStates; ++j) {
			*value += fabs(U->ptr[j] - U2.ptr[j]);
		}
]], {
	eqn = self,
})},
	}
end


-- these are calculated based on cell-centered (or extrapolated) conserved vars
-- they are used to calculate the eigensystem at a cell center or edge 
GLM_MHD.roeVars = table{
	{rho = 'real'},
	{v = 'real3'},
	{hTotal = 'real'},
	{B = 'real3'},
	{X = 'real'},
	{Y = 'real'},
}


-- here's the variables that an eigensystem uses to compute a left, right, or flux transform 
GLM_MHD.eigenVars = table(GLM_MHD.roeVars):append{

	{hHydro = 'real'},
	{aTildeSq = 'real'},

	{Cs = 'real'},
	{CAx = 'real'},
	{Cf = 'real'},
}:append(GLM_MHD.useFixedCh and {} or {
	{Ch = 'real'},
}):append{
	{BStarPerpLen = 'real'},
	{betaY = 'real'},
	{betaZ = 'real'},
	{betaStarY = 'real'},
	{betaStarZ = 'real'},
	{betaStarSq = 'real'},

	{alphaF = 'real'},
	{alphaS = 'real'},

	{sqrtRho = 'real'},
	{sbx = 'real'},
	{Qf = 'real'},
	{Qs = 'real'},
	{Af = 'real'},
	{As = 'real'},
}


function GLM_MHD:getEigenTypeCode()
	return table{
		makestruct.makeStruct('Roe_t', self.roeVars),
		GLM_MHD.super.getEigenTypeCode(self),
	}:concat'\n'
end

function GLM_MHD:eigenWaveCode(side, eig, x, waveIndex)
	eig = '('..eig..')'
	return ({
		eig..'.v.x - '..eig..'.Cf',
		eig..'.v.x - '..eig..'.CAx',
		eig..'.v.x - '..eig..'.Cs',
		eig..'.v.x',
		eig..'.v.x + '..eig..'.Cs',
		eig..'.v.x + '..eig..'.CAx',
		eig..'.v.x + '..eig..'.Cf',
		
		--#warning there's a few PLM routines that expect eigenvalues to be ordered ... so replace them with a eigen_calcMinMaxWaves
		self.useFixedCh and '-Ch' or '-'..eig..'.Ch',
		self.useFixedCh and 'Ch' or eig..'.Ch',
	})[waveIndex+1]
end

-- this is all temporary fix until I properly code the inlining

GLM_MHD.hasCalcDTCode = true

function GLM_MHD:consWaveCodePrefix(side, U, x)
	return template([[
	range_t lambda = calcCellMinMaxEigenvalues_<?=side?>(solver, &<?=U?>, <?=x?>); 
]], {
		side = side,
		U = '('..U..')',
		x = x,
	})
end

function GLM_MHD:consMinWaveCode(side, U, x)
	return 'lambda.min'
end

function GLM_MHD:consMaxWaveCode(side, U, x)
	return 'lambda.max'
end

return GLM_MHD

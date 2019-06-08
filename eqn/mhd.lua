local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local makestruct = require 'eqn.makestruct'
local Equation = require 'eqn.eqn'

local MHD = class(Equation)

MHD.name = 'MHD'

MHD.numWaves = 7
MHD.numIntStates = 8
MHD.numStates = 10

MHD.primVars = table{
	{name='rho', type='real'},
	{name='v', type='real3'},
	{name='P', type='real'},
	{name='B', type='real3'},
	{name='BPot', type='real'},
	{name='ePot', type='real'},	-- for selfgrav
}

MHD.consVars = table{
	{name='rho', type='real'},
	{name='m', type='real3'},
	{name='ETotal', type='real'},
	{name='B', type='real3'},
	{name='BPot', type='real'},
	{name='ePot', type='real'},	-- for selfgrav
}

MHD.mirrorVars = {{'m.x', 'B.x'}, {'m.y', 'B.y'}, {'m.z', 'B.z'}}

MHD.hasEigenCode = true
MHD.hasFluxFromConsCode = true
MHD.roeUseFluxFromCons = true

-- for connections
MHD.useSourceTerm = true

-- hmm, we want init.euler and init.mhd here ...
MHD.initStates = require 'init.euler'

function MHD:init(args)
	MHD.super.init(self, args)

	if self.solver.dim > 1 then
		local NoDiv = require 'op.nodiv'
		self.solver.ops:insert(NoDiv{solver=self.solver})
	end

	-- hmm...
	local SelfGrav = require 'op.selfgrav'
	self.gravOp = SelfGrav{solver=self.solver}
	self.solver.ops:insert(self.gravOp)
end

MHD.guiVars = {
	{name='heatCapacityRatio', value=2},	-- 5/3 for most problems, but 2 for Brio-Wu, so I will just set it here for now (in case something else is reading it before it is set there)
	{name='mu0', value=1},	-- this should be 4 pi for natural units, but I haven't verified that all mu0's are where they should be ...
}

function MHD:getCommonFuncCode()
	return template([[
static inline real calc_eKin(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.v, x); }
static inline real calc_EKin(<?=eqn.prim_t?> W, real3 x) { return W.rho * calc_eKin(W, x); }
static inline real calc_EInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return W.P / (solver->heatCapacityRatio - 1.); }
static inline real calc_eInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_EInt(solver, W) / W.rho; }
static inline real calc_EMag(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.B, x); }
static inline real calc_eMag(<?=eqn.prim_t?> W, real3 x) { return calc_EMag(W, x) / W.rho; }
static inline real calc_PMag(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.B, x); }
static inline real calc_EHydro(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { return calc_EKin(W, x) + calc_EInt(solver, W); }
static inline real calc_eHydro(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { return calc_EHydro(solver, W, x) / W.rho; }
static inline real calc_ETotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { return calc_EKin(W, x) + calc_EInt(solver, W) + calc_EMag(W, x) + W.rho * W.ePot; }
static inline real calc_eTotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { return calc_ETotal(solver, W, x) / W.rho; }
static inline real calc_H(constant <?=solver.solver_t?>* solver, real P) { return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)); }
static inline real calc_h(constant <?=solver.solver_t?>* solver, real rho, real P) { return calc_H(solver, P) / rho; }
static inline real calc_HTotal(<?=eqn.prim_t?> W, real ETotal, real3 x) { return W.P + calc_PMag(W, x) + ETotal; }
static inline real calc_hTotal(<?=eqn.prim_t?> W, real ETotal, real3 x) { return calc_HTotal(W, ETotal, x) / W.rho; }

//notice, this is speed of sound, to match the name convention of eqn/euler
//but Cs in eigen_t is the slow speed
//most the MHD papers use 'a' for the speed of sound
static inline real calc_Cs(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return sqrt(solver->heatCapacityRatio * W.P / W.rho); }
]], {
		solver = self.solver,
		eqn = self,
	})
end

function MHD:getPrimConsCode()
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
	real vSq = coordLenSq(W.v, x);
	real BSq = coordLenSq(W.B, x);
	real EKin = .5 * U.rho * vSq;
	real EMag = .5 * BSq;
	real EPot = U.rho * U.ePot;
	real EInt = U.ETotal - EKin - EMag - EPot;
	W.P = EInt * (solver->heatCapacityRatio - 1.);
	W.P = max(W.P, (real)1e-7);
	W.rho = max(W.rho, (real)1e-7);
	W.BPot = U.BPot;
	W.ePot = U.ePot;
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
	real vSq = coordLenSq(W.v, x);
	real BSq = coordLenSq(W.B, x);
	real EKin = .5 * W.rho * vSq;
	real EMag = .5 * BSq;
	real EInt = W.P / (solver->heatCapacityRatio - 1.);
	U.ETotal = EInt + EKin + EMag;
	U.BPot = W.BPot;
	U.ePot = W.ePot;
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
			+ W.P / (solver->heatCapacityRatio - 1.)
			+ WA.rho * W.ePot,
		.BPot = W.BPot,
		.ePot = W.ePot,
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
			+ U.ETotal
			- WA.rho * U.ePot),
		.BPot = U.BPot,
		.ePot = U.ePot,
	};
}
]], {
		solver = self.solver,
		eqn = self,
	})
end

MHD.initStateCode = [[
<? local xNames = require 'common'().xNames ?>
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
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
	real3 B = real3_zero;

	<?=code?>
	
	<?=eqn.prim_t?> W = {
		.rho = rho,
		.v = cartesianToCoord(v, x),
		.P = P,
		.B = cartesianToCoord(B, x),
		.BPot = 0,
	};
	UBuf[index] = consFromPrim(solver, W, x);
}
]]

MHD.solverCodeFile = 'eqn/mhd.cl'

MHD.displayVarCodeUsesPrims = true

MHD.predefinedDisplayVars = {
	'U rho',
	'U m',
	'U ETotal',
	'U B',
	'U div B',
}

function MHD:getDisplayVars()
	return MHD.super.getDisplayVars(self):append{
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
		{['BPot'] = '*value = U->BPot;'},
		{P = '*value = W.P;'},
		--{PMag = '*value = calc_PMag(W, x);'},
		--{PTotal = '*value = W.P + calc_PMag(W, x);'},
		--{eInt = '*value = calc_eInt(solver, W);'},
		{EInt = '*value = calc_EInt(solver, W);'},
		--{eKin = '*value = calc_eKin(W, x);'},
		{EKin = '*value = calc_EKin(W, x);'},
		--{eHydro = '*value = calc_eHydro(solver, W, x);'},
		{EHydro = '*value = calc_EHydro(solver, W, x);'},
		--{eMag = '*value = calc_eMag(W, x);'},
		{EMag = '*value = calc_EMag(W, x);'},
		--{eTotal = '*value = U->ETotal / W.rho;'},
		{S = '*value = W.P / pow(W.rho, (real)solver->heatCapacityRatio);'},
		{H = '*value = calc_H(solver, W.P);'},
		--{h = '*value = calc_H(solver, W.P) / W.rho;'},
		--{HTotal = '*value = calc_HTotal(W, U->ETotal, x);'},
		--{hTotal = '*value = calc_hTotal(W, U->ETotal, x);'},
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
MHD.roeVars = table{
	{name='rho', type='real'},
	{name='v', type='real3'},
	{name='hTotal', type='real'},
	{name='B', type='real3'},
	{name='X', type='real'},
	{name='Y', type='real'},
}


-- here's the variables that an eigensystem uses to compute a left, right, or flux transform 
MHD.eigenVars = table(MHD.roeVars):append{

	{name='hHydro', type='real'},
	{name='aTildeSq', type='real'},

	{name='Cs', type='real'},
	{name='CAx', type='real'},
	{name='Cf', type='real'},

	{name='BStarPerpLen', type='real'},
	{name='betaY', type='real'},
	{name='betaZ', type='real'},
	{name='betaStarY', type='real'},
	{name='betaStarZ', type='real'},
	{name='betaStarSq', type='real'},

	{name='alphaF', type='real'},
	{name='alphaS', type='real'},

	{name='sqrtRho', type='real'},
	{name='sbx', type='real'},
	{name='Qf', type='real'},
	{name='Qs', type='real'},
	{name='Af', type='real'},
	{name='As', type='real'},
}


function MHD:getEigenTypeCode()
	return table{
		makestruct.makeStruct('Roe_t', self.roeVars),
		MHD.super.getEigenTypeCode(self),
	}:concat'\n'
end

function MHD:eigenWaveCode(side, eig, x, waveIndex)
	eig = '('..eig..')'
	return ({
		eig..'.v.x - '..eig..'.Cf',
		eig..'.v.x - '..eig..'.CAx',
		eig..'.v.x - '..eig..'.Cs',
		eig..'.v.x',
		eig..'.v.x + '..eig..'.Cs',
		eig..'.v.x + '..eig..'.CAx',
		eig..'.v.x + '..eig..'.Cf',
	})[waveIndex+1]
end

-- this is all temporary fix until I properly code the inlining

MHD.hasCalcDTCode = true

function MHD:consWaveCodePrefix(side, U, x)
	return template([[
	range_t lambda = calcCellMinMaxEigenvalues_<?=side?>(solver, &<?=U?>, <?=x?>); 
]], {
		side = side,
		U = '('..U..')',
		x = x,
	})
end

function MHD:consMinWaveCode(side, U, x)
	return 'lambda.min'
end

function MHD:consMaxWaveCode(side, U, x)
	return 'lambda.max'
end

return MHD

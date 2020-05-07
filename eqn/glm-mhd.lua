local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local template = require 'template'
local makestruct = require 'eqn.makestruct'
local constants = require 'constants'
local Equation = require 'eqn.eqn'

local MHD = class(Equation)

MHD.name = 'GLM-MHD'

MHD.numWaves = 9
MHD.numIntStates = 9
MHD.numStates = 10

MHD.primVars = table{
	{name='rho', type='real', units='kg/m^3'},
	{name='v', type='real3', units='m/s', variance='u'},
	{name='P', type='real', units='kg/(m*s^2)'},
	{name='B', type='real3', units='kg/(C*s)', variance='l'},
	{name='psi', type='real', units='kg/(C*s)'},
	{name='ePot', type='real', units='m^2/s^2'},
}

MHD.consVars = table{
	{name='rho', type='real', units='kg/m^3'},
	{name='m', type='real3', units='kg/(m^2*s)', variance='u'},
	{name='ETotal', type='real', units='kg/(m*s^2)'},
	{name='B', type='real3', units='kg/(C*s)', variance='l'},
	{name='psi', type='real', units='kg/(C*s)'},
	{name='ePot', type='real', units='m^2/s^2'},
}

MHD.hasEigenCode = true
MHD.hasFluxFromConsCode = true
MHD.roeUseFluxFromCons = true
MHD.useSourceTerm = true
MHD.useConstrainU = true

-- TODO this is broken
MHD.useFixedCh = true -- true = use a gui var, false = calculate by max(|v_i|+Cf)

-- hmm, we want init.euler and init.mhd here ...
MHD.initStates = require 'init.euler'

function MHD:init(args)
	MHD.super.init(self, args)

	local UpdatePsi = require 'op.glm-mhd-update-psi'
	self.solver.ops:insert(UpdatePsi{solver=self.solver})
	
	local SelfGrav = require 'op.selfgrav'
	self.gravOp = SelfGrav{solver=self.solver}
	self.solver.ops:insert(self.gravOp)
end

--[[
B^2 = kg^2/(C^2 s^2)
B^2 / mu = kg^2/(C^2 s^2) * C^2/(kg*m) = kg/(m s^2)
--]]
MHD.guiVars = table{
	{name='heatCapacityRatio', value=2},	-- 5/3 for most problems, but 2 for Brio-Wu, so I will just set it here for now (in case something else is reading it before it is set there)
	
	-- works good as mu0 = 1
	-- solver->mu0 / unit_kg_m_per_C2 == 1
	-- solver->mu0 / (unit_kg * unit_m / (unit_C * unit_C)) == 1
	-- unit_C = sqrt((unit_kg * unit_m) / solver->mu0);
	{name='mu0', value=constants.vacuumPermeability_in_kg_m_per_C2, units='(kg*m)/C^2'},
	
	{name='Cp', value=1, compileTime=true},
}

if MHD.useFixedCh then
	MHD.guiVars:insert{name='Ch', value=.1}
end

function MHD:getCommonFuncCode()
	return template([[
real calc_eKin(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.v, x); }
real calc_EKin(<?=eqn.prim_t?> W, real3 x) { return W.rho * calc_eKin(W, x); }
real calc_EInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return W.P / (solver->heatCapacityRatio - 1.); }
real calc_eInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_EInt(solver, W) / W.rho; }
real calc_EM_energy(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.B, x) / (solver->mu0 / unit_kg_m_per_C2); }
real calc_PMag(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.B, x) / (solver->mu0 / unit_kg_m_per_C2); }
real calc_EHydro(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { return calc_EKin(W, x) + calc_EInt(solver, W); }
real calc_eHydro(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { return calc_EHydro(solver, W, x) / W.rho; }
real calc_ETotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { return calc_EKin(W, x) + calc_EInt(solver, W) + calc_EM_energy(solver, W, x); }
real calc_eTotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { return calc_ETotal(solver, W, x) / W.rho; }
real calc_H(constant <?=solver.solver_t?>* solver, real P) { return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)); }
real calc_h(constant <?=solver.solver_t?>* solver, real rho, real P) { return calc_H(solver, P) / rho; }
real calc_HTotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real ETotal, real3 x) { return W.P + calc_PMag(solver, W, x) + ETotal; }
real calc_hTotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real ETotal, real3 x) { return calc_HTotal(solver, W, ETotal, x) / W.rho; }

//notice, this is speed of sound, to match the name convention of eqn/euler
//but Cs in eigen_t is the slow speed
//most the MHD papers use 'a' for the speed of sound
real calc_Cs(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return sqrt(solver->heatCapacityRatio * W.P / W.rho); }

]], {
		solver = self.solver,
		eqn = self,
	})
end

function MHD:getPrimConsCode()
	return template([[
<?=eqn.prim_t?> primFromCons(
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
	real EMag = .5 * BSq / (solver->mu0 / unit_kg_m_per_C2);
	real EInt = U.ETotal - EKin - EMag;
	W.P = EInt * (solver->heatCapacityRatio - 1.);
	W.P = max(W.P, (real)1e-7);
	W.rho = max(W.rho, (real)1e-7);
	W.psi = U.psi;
	W.ePot = U.ePot;
	return W;
}

<?=eqn.cons_t?> consFromPrim(
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
	real EMag = .5 * BSq / (solver->mu0 / unit_kg_m_per_C2);
	real EInt = W.P / (solver->heatCapacityRatio - 1.);
	U.ETotal = EInt + EKin + EMag;
	U.psi = W.psi;
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
			+ real3_dot(W.B, WA.B) / (solver->mu0 / unit_kg_m_per_C2)
			+ W.P / (solver->heatCapacityRatio - 1.),
		.psi = W.psi,
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
			- real3_dot(U.B, WA.B) / (solver->mu0 / unit_kg_m_per_C2)
			+ U.ETotal),
		.psi = U.psi,
		.ePot = U.ePot,
	};
}
]], {
		solver = self.solver,
		eqn = self,
	})
end

MHD.initStateCode = [[
<? local xNames = require 'common'.xNames ?>
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
	real ePot = 0;

	<?=code?>
	
	<?=eqn.prim_t?> W = {
		.rho = rho,
		.v = cartesianToCoord(v, x),
		.P = P,
		.B = cartesianToCoord(B, x),
		.psi = 0,
		.ePot = ePot,
	};
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

function MHD:getInitStateCode()

	-- where do I put this to make it the default value for MHD solvers,
	-- but not override a value set by the init state?
	self.guiVars.coulomb.value = math.sqrt(self.guiVars.kilogram.value * self.guiVars.meter.value / self.guiVars.mu0.value)

	return MHD.super.getInitStateCode(self)
end

MHD.solverCodeFile = 'eqn/glm-mhd.cl'

MHD.displayVarCodeUsesPrims = true

MHD.predefinedDisplayVars = {
	'U rho',
	'U m mag',
	'U ETotal',
	'U B mag',
	'U div B',
	'U psi',
}

function MHD:getDisplayVars()
	return MHD.super.getDisplayVars(self):append{
		{name='v', code='value.vreal3 = W.v;', type='real3', units='m/s'},
		{name='div B', code=template([[
	value.vreal = .5 * (0.
<? 
for j=0,solver.dim-1 do 
?>		+ (U[solver->stepsize.s<?=j?>].<?=field?>.s<?=j?> 
			- U[-solver->stepsize.s<?=j?>].<?=field?>.s<?=j?>
		) / solver->grid_dx.s<?=j?>
<? 
end 
?>	);
]], {solver=self.solver, field='B'}), units='kg/(C*m*s)'},
		{name='P', code='value.vreal = W.P;', units='kg/(m*s^2)'},
		--{name='PMag', code='value.vreal = calc_PMag(solver, W, x);'},
		--{name='PTotal', code='value.vreal = W.P + calc_PMag(solver, W, x);'},
		--{name='eInt', code='value.vreal = calc_eInt(solver, W);'},
		{name='EInt', code='value.vreal = calc_EInt(solver, W);'},
		--{name='eKin', code='value.vreal = calc_eKin(W, x);'},
		{name='EKin', code='value.vreal = calc_EKin(W, x);'},
		--{name='eHydro', code='value.vreal = calc_eHydro(solver, W, x);'},
		{name='EHydro', code='value.vreal = calc_EHydro(solver, W, x);'},
		{name='EM energy', code='value.vreal = calc_EM_energy(solver, W, x);'},
		--{name='eTotal', code='value.vreal = U->ETotal / W.rho;'},
		{name='S', code='value.vreal = W.P / pow(W.rho, (real)solver->heatCapacityRatio);'},
		{name='H', code='value.vreal = calc_H(solver, W.P);'},
		--{name='h', code='value.vreal = calc_H(solver, W.P) / W.rho;'},
		--{name='HTotal', code='value.vreal = calc_HTotal(solver, W, U->ETotal, x);'},
		--{name='hTotal', code='value.vreal = calc_hTotal(solver, W, U->ETotal, x);'},
		--{name='Cs', code='value.vreal = calc_Cs(solver, W); },
		{name='primitive reconstruction error', code=template([[
		//prim have just been reconstructed from cons
		//so reconstruct cons from prims again and calculate the difference
		<?=eqn.cons_t?> U2 = consFromPrim(solver, W, x);
		value.vreal = 0;
		for (int j = 0; j < numIntStates; ++j) {
			value.vreal += fabs(U->ptr[j] - U2.ptr[j]);
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
	{name='Ch', type='real'},
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

function MHD:eigenWaveCode(n, eig, x, waveIndex)
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
		'-'..eig..'.Ch',
		eig..'.Ch',
	})[waveIndex+1] or error("got a bad waveIndex")
end

function MHD:consWaveCodePrefix(n, U, x)
	return template([[
	range_t lambda = calcCellMinMaxEigenvalues(solver, &<?=U?>, <?=x?>, <?=n?>); 
]], {
		n = n,
		U = '('..U..')',
		x = x,
	})
end

function MHD:consMinWaveCode(n, U, x)
	return 'lambda.min'
end

function MHD:consMaxWaveCode(n, U, x)
	return 'lambda.max'
end

return MHD

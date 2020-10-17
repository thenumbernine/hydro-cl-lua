local class = require 'ext.class'
local table = require 'ext.table'
local constants = require 'hydro.constants'
local Struct = require 'hydro.code.struct'
local Equation = require 'hydro.eqn.eqn'

local MHD = class(Equation)

MHD.name = 'MHD'

MHD.numWaves = 7
MHD.numIntStates = 8

MHD.roeUseFluxFromCons = true

-- for connections
MHD.useSourceTerm = true
MHD.useConstrainU = true

-- hmm, we want init.euler and init.mhd here ...
MHD.initConds = require 'hydro.init.euler':getList()


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

function MHD:init(args)
	
	-- TODO redo the mhd equations for a background grid metric, and take note of covariance/contravariance
	self.primVars = table{
		{name='rho', type='real', units='kg/m^3'},
		{name='v', type='real3', units='m/s', variance='u'},
		{name='P', type='real', units='kg/(m*s^2)'},
		{name='B', type='real3', units='kg/(C*s)', variance='l'},
		{name='psi', type='real', units='kg/(C*s)'},
		{name='ePot', type='real', units='m^2/s^2'},
	}

	self.consVars = table{
		{name='rho', type='real', units='kg/m^3'},
		{name='m', type='real3', units='kg/(m^2*s)', variance='u'},
		{name='ETotal', type='real', units='kg/(m*s^2)'},
		{name='B', type='real3', units='kg/(C*s)', variance='l'},
		{name='psi', type='real', units='kg/(C*s)'},
		{name='ePot', type='real', units='m^2/s^2'},
	}
	
	if args.incompressible then
		self.consVars:insert{name='mPot', type='real', units='kg/(m*s)'}
		self.primVars:insert{name='mPot', type='real', units='kg/(m*s)'}
	end

	MHD.super.init(self, args)
	
	local solver = self.solver

	self.roeStruct = Struct{solver=solver, name='roe_t', vars=self.roeVars}
	self.roeStruct:makeType()
	self.roe_t = self.roeStruct.typename

	if require 'hydro.solver.meshsolver'.is(solver) then
		print("not using ops (selfgrav, nodiv, etc) with mesh solvers yet")
	else
		if solver.dim > 1 then
			local NoDiv = require 'hydro.op.nodiv'()
			solver.ops:insert(NoDiv{
				solver = solver,
				potentialField = 'psi',
			})
		end

		local SelfGrav = require 'hydro.op.selfgrav'
		self.gravOp = SelfGrav{solver=solver}
		solver.ops:insert(self.gravOp)
	
		if args.incompressible then
			local NoDiv = require 'hydro.op.nodiv'{
				poissonSolver = require 'hydro.op.poisson_jacobi',	-- krylov is having errors.  TODO bug in its boundary code?
			}
			self.solver.ops:insert(NoDiv{
				solver = self.solver,
				vectorField = 'm',
				potentialField = 'mPot',
			
				-- div v = 0
				-- div (m/ρ) = 0
				-- 1/ρ div m - 1/ρ^2 m dot grad ρ = 0
				-- div m = (m dot grad ρ)/ρ 
				chargeCode = self:template[[
	<? for j=0,solver.dim-1 do ?>{
		global const <?=eqn.cons_t?>* Ujm = U - solver->stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + solver->stepsize.s<?=j?>;
		real drho_dx = (Ujp->rho - Ujm->rho) * (.5 / solver->grid_dx.s<?=j?>);
		source -= drho_dx * U->m.s<?=j?> / U->rho;
	}<? end ?>
]],
			})
		end
	end
end

--[[
B^2 = kg^2/(C^2 s^2)
B^2 / mu = kg^2/(C^2 s^2) * C^2/(kg*m) = kg/(m s^2)
--]]
MHD.guiVars = {
	{name='heatCapacityRatio', value=2},	-- 5/3 for most problems, but 2 for Brio-Wu, so I will just set it here for now (in case something else is reading it before it is set there)
	
	-- works good as mu0 = 1
	-- solver->mu0 / unit_kg_m_per_C2 == 1
	-- solver->mu0 / (unit_kg * unit_m / (unit_C * unit_C)) == 1
	-- unit_C = sqrt((unit_kg * unit_m) / solver->mu0);
	{name='mu0', value=constants.vacuumPermeability_in_kg_m_per_C2, units='(kg*m)/C^2'},
}

function MHD:initCodeModules()
	MHD.super.initCodeModules(self)
	
	-- TODO find a better place to put this
	
	local solver = self.solver

	solver.modules:add{
		name = 'roe_t',
		structs = {self.roeStruct},
	}
	solver.solverModulesEnabled.roe_t = true

	-- added by request only, so I don't have to compile the real3x3 code
	solver.modules:add{
		name = 'calcCellMinMaxEigenvalues',
		code = self:template[[
//called from calcDT
range_t calcCellMinMaxEigenvalues(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	<?=eqn.cons_t?> U_ = cons_rotateFrom(U, n);
	<?=eqn.prim_t?> W = primFromCons(solver, U_, x);
	
#if 0
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real3 v = W.v;
	real3 B = W.B;
	
	real BSq = coordLenSq(B, x);
	real invRho = 1./W.rho;
	
	real aSq = solver->heatCapacityRatio * W.P * invRho;
	real B_n = normal_vecDotN1(n, B);
	real CaxSq = B_n * B_n * invRho;
	real CaSq = BSq * invRho;
	
	real CStarSq = .5 * (CaSq + aSq);
	real sqrtCfsDiscr = sqrt(max(0., CStarSq * CStarSq - aSq * CaxSq));
	
	real CfSq = CStarSq + sqrtCfsDiscr;
	real CsSq = CStarSq - sqrtCfsDiscr;

	real Cf = sqrt(CfSq);
	real Cs = sqrt(max(CsSq, 0.));
	real v_n = normal_vecDotN1(n, v);
	return (range_t){.min=v_n - Cf, .max=v_n + Cf};
#else
	const real gamma = solver->heatCapacityRatio;
	const real gamma_1 = gamma - 1.;
	const real gamma_2 = gamma - 2.;
	
	real rho = W.rho;
	real3 v = W.v;
	real3 B = W.B;
	real hTotal = .5 * coordLenSq(W.v, x) + (W.P * gamma / gamma_1 + coordLenSq(B, x)) / W.rho;

	//the rest of this matches calcEigenBasis:

	real _1_rho = 1. / rho;
	real vSq = coordLenSq(v, x);
#warning consider g_ij	
	real BPerpSq = B.y*B.y + B.z*B.z;
	real BStarPerpSq = (gamma_1 - gamma_2) * BPerpSq;
	real CAxSq = B.x*B.x*_1_rho;
	real CASq = CAxSq + BPerpSq * _1_rho;
	real hHydro = hTotal - CASq;
	// hTotal = (EHydro + EMag + P)/rho
	// hHydro = hTotal - CASq, CASq = EMag/rho
	// hHydro = eHydro + P/rho = eKin + eInt + P/rho
	// hHydro - eKin = eInt + P/rho = (1./(gamma-1) + 1) P/rho = gamma/(gamma-1) P/rho
	// a^2 = (gamma-1)(hHydro - eKin) = gamma P / rho
	real aTildeSq = max((gamma_1 * (hHydro - .5 * vSq) - gamma_2), 1e-20);

	real BStarPerpSq_rho = BStarPerpSq * _1_rho;
	real CATildeSq = CAxSq + BStarPerpSq_rho;
	real CStarSq = .5 * (CATildeSq + aTildeSq);
	real CA_a_TildeSqDiff = .5 * (CATildeSq - aTildeSq);
	real sqrtDiscr = sqrt(CA_a_TildeSqDiff * CA_a_TildeSqDiff + aTildeSq * BStarPerpSq_rho);
	
	real CfSq = CStarSq + sqrtDiscr;
	real Cf = sqrt(CfSq);

	real CsSq = aTildeSq * CAxSq / CfSq;
	real Cs = sqrt(CsSq);

	real lambdaFastMin = v.x - Cf;
	real lambdaFastMax = v.x + Cf;
	
	return (range_t){
		.min = lambdaFastMin,
		.max = lambdaFastMax,
	};
#endif
}
]],
	}

	-- TODO don't put this here, instead make it a depends of the calcDT/consWaveCodePrefix code below that references it.
	solver.solverModulesEnabled.calcCellMinMaxEigenvalues = true
end

function MHD:initCodeModule_fluxFromCons()
	self.solver.modules:add{
		name = 'fluxFromCons',
		depends = {
			'solver.solver_t',
			'eqn.cons_t',
			'eqn.prim_t',
			'eqn.prim-cons',	-- primFromCons
			'normal_t',
			'coordLenSq',
		},
		code = self:template[[
<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	<?=eqn.prim_t?> W = primFromCons(solver, U, x);
	real vj = normal_vecDotN1(n, W.v);
	real Bj = normal_vecDotN1(n, W.B);
	real BSq = coordLenSq(W.B, x);
	real BDotV = real3_dot(W.B, W.v);
	real PMag = .5 * BSq / (solver->mu0 / unit_kg_m_per_C2);
	real PTotal = W.P + PMag;
	real HTotal = U.ETotal + PTotal;
	
	<?=eqn.cons_t?> F;
	F.rho = normal_vecDotN1(n, U.m);
	F.m = real3_sub(real3_real_mul(U.m, vj), real3_real_mul(U.B, Bj / (solver->mu0 / unit_kg_m_per_C2)));
	F.m.x += PTotal * normal_u1x(n);
	F.m.y += PTotal * normal_u1y(n);
	F.m.z += PTotal * normal_u1z(n);
	F.B = real3_sub(real3_real_mul(U.B, vj), real3_real_mul(W.v, Bj));
	F.ETotal = HTotal * vj - BDotV * Bj / (solver->mu0 / unit_kg_m_per_C2);
	F.psi = 0;
	F.ePot = 0;
	return F;
}
]],
	}
end

function MHD:initCodeModuleCommon()
	self.solver.modules:add{
		name = 'eqn.common',
		code = self:template[[
real calc_eKin(<?=eqn.prim_t?> W, real3 x) { 
	return .5 * coordLenSq(W.v, x);
}

real calc_EKin(<?=eqn.prim_t?> W, real3 x) { 
	return W.rho * calc_eKin(W, x); 
}

real calc_EInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { 
	return W.P / (solver->heatCapacityRatio - 1.); 
}

real calc_eInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { 
	return calc_EInt(solver, W) / W.rho; 
}

//units: 
//B has units kg/(C*s)
//mu0 has units kg*m/C^2
//PMag = 1/2 B^2 / mu0 has units kg/(m*s^2)
real calc_EM_energy(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { 
	return .5 * coordLenSq(W.B, x) / (solver->mu0 / unit_kg_m_per_C2); 
}

//same as calc_EM_energy
real calc_PMag(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { 
	return .5 * coordLenSq(W.B, x) / (solver->mu0 / unit_kg_m_per_C2); 
}

real calc_EHydro(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { 
	return calc_EKin(W, x) + calc_EInt(solver, W); 
}

real calc_eHydro(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { 
	return calc_EHydro(solver, W, x) / W.rho; 
}

real calc_ETotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { 
	return calc_EKin(W, x) + calc_EInt(solver, W) + calc_EM_energy(solver, W, x); 
}

real calc_eTotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) { 
	return calc_ETotal(solver, W, x) / W.rho; 
}

real calc_H(constant <?=solver.solver_t?>* solver, real P) { 
	return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)); 
}

real calc_h(constant <?=solver.solver_t?>* solver, real rho, real P) { 
	return calc_H(solver, P) / rho; 
}

real calc_HTotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real ETotal, real3 x) { 
	return W.P + calc_PMag(solver, W, x) + ETotal; 
}

real calc_hTotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real ETotal, real3 x) { 
	return calc_HTotal(solver, W, ETotal, x) / W.rho; 
}

//notice, this is speed of sound, to match the name convention of hydro/eqn/euler
//but Cs in eigen_t is the slow speed
//most the MHD papers use 'a' for the speed of sound
real calc_Cs(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { 
	return sqrt(solver->heatCapacityRatio * W.P / W.rho);
}

//CA = B/sqrt(mu0 rho)
//B has units kg/(C*s)
//mu0 has units kg*m/C^2
//rho has units kg/m^3
//CA has units m/s
real3 calc_CA(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U) {
	return real3_real_mul(U.B, 1./sqrt(U.rho * solver->mu0 / unit_kg_m_per_C2));
}
]],
	}
end


function MHD:getModuleDependsSolver() 
	return {'eqn.prim-cons'}
end

function MHD:initCodeModulePrimCons()
	self.solver.modules:add{
		name = 'eqn.prim-cons',
		depends = {
			'real3',
			'solver.solver_t',
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
]],
	}
end

MHD.initCondCode = [[
kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	
	global <?=eqn.cons_t?>* U = UBuf + index;
	
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
	//ignored:
	real3 D = real3_zero;

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
]]

function MHD:getInitCondCode()

	-- where do I put this to make it the default value for MHD solvers,
	-- but not override a value set by the init state?
	self.guiVars.coulomb.value = math.sqrt(self.guiVars.kilogram.value * self.guiVars.meter.value / self.guiVars.mu0.value)

	return MHD.super.getInitCondCode(self)
end

MHD.solverCodeFile = 'hydro/eqn/mhd.cl'

MHD.displayVarCodeUsesPrims = true

MHD.predefinedDisplayVars = {
	'U rho',
	'U m',
	'U ETotal',
	'U D',
	'U B',
	'U div B',
}

function MHD:getDisplayVars()
	local vars = MHD.super.getDisplayVars(self)
	vars:append{
		{name='v', code='value.vreal3 = W.v;', type='real3', units='m/s'},
		{name='P', code='value.vreal = W.P;', units='kg/(m*s^2)'},
		{name='PMag', code='value.vreal = calc_PMag(solver, W, x);', units='kg/(m*s^2)'},
		{name='PTotal', code='value.vreal = W.P + calc_PMag(solver, W, x);', units='kg/(m*s^2)'},
		{name='eInt', code='value.vreal = calc_eInt(solver, W);', units='m^2/s^2'},
		{name='EInt', code='value.vreal = calc_EInt(solver, W);', units='kg/(m*s^2)'},
		{name='eKin', code='value.vreal = calc_eKin(W, x);', units='m^2/s^2'},
		{name='EKin', code='value.vreal = calc_EKin(W, x);', units='kg/(m*s^2)'},
		{name='eHydro', code='value.vreal = calc_eHydro(solver, W, x);', units='m^2/s^2'},
		{name='EHydro', code='value.vreal = calc_EHydro(solver, W, x);', units='kg/(m*s^2)'},
		{name='EM energy', code='value.vreal = calc_EM_energy(solver, W, x);', units='kg/(m*s^2)'},
		{name='eTotal', code='value.vreal = U->ETotal / W.rho;', units='m^2/s^2'},
		{name='S', code='value.vreal = W.P / pow(W.rho, (real)solver->heatCapacityRatio);'},
		{name='H', code='value.vreal = calc_H(solver, W.P);', units='kg/(m*s^2)'},
		{name='h', code='value.vreal = calc_H(solver, W.P) / W.rho;', units='m^2/s^2'},
		{name='HTotal', code='value.vreal = calc_HTotal(solver, W, U->ETotal, x);', units='kg/(m*s^2)'},
		{name='hTotal', code='value.vreal = calc_hTotal(solver, W, U->ETotal, x);', units='m^2/s^2'},
		{name='speed of sound', code='value.vreal = calc_Cs(solver, W);', units='m/s'},
		{name='alfven velocity', code='value.vreal3 = calc_CA(solver, *U);', type='real3', units='m/s'},
		{name='Mach number', code='value.vreal = coordLen(W.v, x) / calc_Cs(solver, W);'},
		{name='temperature', code=self:template[[
<? local clnumber = require 'cl.obj.number' ?>
<? local materials = require 'hydro.materials' ?>
#define C_v				<?=('%.50f'):format(materials.Air.C_v)?>
	value.vreal = calc_eInt(solver, W) / C_v;
]], units='K'},
		{name='primitive reconstruction error', code=self:template[[
	//prim have just been reconstructed from cons
	//so reconstruct cons from prims again and calculate the difference
	<?=eqn.cons_t?> U2 = consFromPrim(solver, W, x);
	value.vreal = 0;
	for (int j = 0; j < numIntStates; ++j) {
		value.vreal += fabs(U->ptr[j] - U2.ptr[j]);
	}
]]},
	}

	if self.gravOp then
		vars:insert{
			name='gravity', 
			code=self:template[[
	if (!OOB(1,1)) {
		value.vreal3 = calcGravityAccel<?=eqn.gravOp.name?>(solver, U);
	}
]], 
			type='real3', 
			units='m/s^2',
		}
	end

	vars:insert(self:createDivDisplayVar{
		field = 'v', 
		getField = function(U, j)
			return U..'->m.s'..j..' / '..U..'->rho'
		end,
		units = '1/s',
	} or nil)

	vars:insert(self:createCurlDisplayVar{
		field = 'v',
		getField = function(U, j)
			return U..'->m.s'..j..' / '..U..'->rho'
		end,
		units = '1/s',
	} or nil)

	vars:insert(self:createDivDisplayVar{field='B', units='kg/(C*m*s)'} or nil)
	vars:insert(self:createCurlDisplayVar{field='B', units='kg/(C*m*s)'} or nil)

	return vars
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
	})[waveIndex+1] or error("got a bad waveIndex")
end

function MHD:consWaveCodePrefix(n, U, x)
	return self:template([[
	range_t lambda = calcCellMinMaxEigenvalues(solver, <?=U?>, <?=x?>, <?=n?>); 
]], {
		n = n,
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

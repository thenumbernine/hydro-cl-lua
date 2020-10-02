local class = require 'ext.class'
local table = require 'ext.table'
local constants = require 'hydro.constants'
local Struct = require 'hydro.code.struct'
local Equation = require 'hydro.eqn.eqn'

local MHD = class(Equation)

MHD.name = 'GLM-MHD'

MHD.numWaves = 9
MHD.numIntStates = 9

MHD.roeUseFluxFromCons = true
MHD.useSourceTerm = true
MHD.useConstrainU = true

-- TODO this is broken
MHD.useFixedCh = true -- true = use a gui var, false = calculate by max(|v_i|+Cf)


MHD.use2009MignoneEqn3 = false	-- 2003 Mignone eqn 3
MHD.use2002DednerEqn24 = true	-- 2002 Dedner eqn 24, EGLM
MHD.use2002DednerEqn38 = true	-- 2002 Dedner eqn 38

-- hmm, we want init.euler and init.mhd here ...
MHD.initConds = require 'hydro.init.euler'


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


function MHD:init(args)
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

	local solver = assert(args.solver)
	if require 'hydro.solver.meshsolver'.is(solver) then
		print("not divergence with mesh solvers yet")
		-- these don't work with mesh solver
		self.use2002DednerEqn24 = false
		self.use2002DednerEqn38 = false
	end

	MHD.super.init(self, args)
	
	self.roeStruct = Struct{solver=solver, name='roe_t', vars=self.roeVars}
	self.roeStruct:makeType()
	self.roe_t = self.roeStruct.typename

	local UpdatePsi = require 'hydro.op.glm-mhd-update-psi'
	solver.ops:insert(UpdatePsi{solver=solver})
	
	if require 'hydro.solver.meshsolver'.is(solver) then
		print("not using ops (selfgrav, nodiv, etc) with mesh solvers yet")
	else
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

function MHD:initCodeModules()
	MHD.super.initCodeModules(self)
	local solver = self.solver

	solver.modules:add{
		name = 'roe_t',
		structs = {self.roeStruct},
	}
	solver.solverModulesEnabled.roe_t = true
end

function MHD:initCodeModule_fluxFromCons()
	self.solver.modules:add{
		name = 'eigen_forCell',
		depends = {
			'solver.solver_t',
			'roe_t',
			'eqn.prim-cons',	-- primFromCons
			'normal_t',
			'coordLenSq',
		},
		code = self:template[[
<?=eqn.eigen_t?> eigen_forCell(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	<?=eqn.prim_t?> W = primFromCons(solver, U, x);
	real PMag = .5 * coordLenSq(W.B, x);
	real hTotal = (U.ETotal + W.P + PMag) / W.rho;
	<?=eqn.roe_t?> roe = {
		.rho = W.rho,
		.v = W.v,
		.hTotal = hTotal,
		.B = W.B,
		.X = 0,
		.Y = 1,
	};
	return eigen_forRoeAvgs(solver, roe, x);
}
]],
	}
	self.solver.solverModulesEnabled['eigen_forCell'] = true

	self.solver.modules:add{
		name = 'fluxFromCons',
		depends = {
			'solver.solver_t',
			'eigen_forCell',
			'normal_t',
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
	
<? if not eqn.useFixedCh then ?>
	//TODO don't need the whole eigen here, just the Ch
	real Ch = 0;
	<? for side=0,solver.dim-1 do ?>{
		<?=eqn.eigen_t?> eig = eigen_forCell(solver, U, x, normal_fromSide<?=side?>(x));
		Ch = max(Ch, eig.Ch);
	}<? end ?>
<? else ?>
	real Ch = solver->Ch;
<? end ?>

	<?=eqn.cons_t?> F;
	F.rho = normal_vecDotN1(n, U.m);
	F.m = real3_sub(real3_real_mul(U.m, vj), real3_real_mul(U.B, Bj / (solver->mu0 / unit_kg_m_per_C2)));
	F.m.x += PTotal * normal_l1x(n);
	F.m.y += PTotal * normal_l1y(n);
	F.m.z += PTotal * normal_l1z(n);
	F.B = real3_sub(real3_real_mul(U.B, vj), real3_real_mul(W.v, Bj));
	F.psi = Ch * Ch;
	F.B.x += F.psi * normal_l1x(n);
	F.B.y += F.psi * normal_l1y(n);
	F.B.z += F.psi * normal_l1z(n);
	F.ETotal = HTotal * vj - BDotV * Bj / (solver->mu0 / unit_kg_m_per_C2);
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

MHD.needsInitDerivs = true
MHD.initCondCode = [[
kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
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

kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cellBuf[index].pos;
	global <?=eqn.cons_t?>* U = UBuf + index;
<? if require 'hydro.solver.meshsolver'.is(solver) then
?>	U->psi = 0.;
<? else
?>	U->psi = .5 * (0.
<? 	for j=0,solver.dim-1 do 
?>		+ (U[solver->stepsize.s<?=j?>].B.s<?=j?> 
			- U[-solver->stepsize.s<?=j?>].B.s<?=j?>
		) / solver->grid_dx.s<?=j?>
<? 	end 
?>	);
<? end
?>
}
]]

function MHD:getInitCondCode()
	-- where do I put this to make it the default value for MHD solvers,
	-- but not override a value set by the init state?
	self.guiVars.coulomb.value = math.sqrt(self.guiVars.kilogram.value * self.guiVars.meter.value / self.guiVars.mu0.value)

	return MHD.super.getInitCondCode(self)
end

MHD.solverCodeFile = 'hydro/eqn/glm-mhd.cl'

MHD.displayVarCodeUsesPrims = true

MHD.predefinedDisplayVars = {
	'U rho',
	'U m',
	'U ETotal',
	'U B',
	'U div B',
	'U psi',
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
		
		--#warning there's a few PLM routines that expect eigenvalues to be ordered ... so replace them with a eigen_calcMinMaxWaves
		'-'..eig..'.Ch',
		eig..'.Ch',
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

function MHD:consMinWaveCode(n, U, x)
	return 'lambda.min'
end

function MHD:consMaxWaveCode(n, U, x)
	return 'lambda.max'
end

return MHD

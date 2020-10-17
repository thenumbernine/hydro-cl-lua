local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'


local Euler = class(Equation)
Euler.name = 'Euler'

-- ePot is the 6th param
-- which means it's now in the derivBuf, but it is always zero
-- so TODO a new variable for deriv size vs cons_t size?
--Euler.numStates = 6	

Euler.numWaves = 5
Euler.numIntStates = 5	-- don't bother integrate ePot

Euler.roeUseFluxFromCons = true

-- the only source term that the Euler equations has is the connection coefficients of the velocity vector
-- maybe later I will automatically flag what elements are vectors
-- and automatically add connection coefficients
--
-- I'm also using the source term for the viscosity ... if you choose to use FANS
Euler.useSourceTerm = true
Euler.useConstrainU = true

Euler.initConds = require 'hydro.init.euler':getList()


function Euler:init(args)
	
	-- TODO primVars doesn't autogen displayVars, and therefore units doesn't matter
	
	self.primVars = table{
		{name='rho', type='real', units='kg/m^3'},
		{name='v', type='real3', units='m/s', variance='u'},			-- contravariant
		{name='P', type='real', units='kg/(m*s^2)'},
		{name='ePot', type='real', units='m^2/s^2'},
	}

	self.consVars = table{
		{name='rho', type='real', units='kg/m^3'},
		{name='m', type='real3', units='kg/(m^2*s)', variance='u'},	-- contravariant
		{name='ETotal', type='real', units='kg/(m*s^2)'},
		{name='ePot', type='real', units='m^2/s^2'},
	}

	if args.incompressible then
		self.consVars:insert{name='mPot', type='real', units='kg/(m*s)'}
		self.primVars:insert{name='mPot', type='real', units='kg/(m*s)'}
	end

	Euler.super.init(self, args)

	if require 'hydro.solver.meshsolver'.is(self.solver) then
		print("not using ops (selfgrav, nodiv, etc) with mesh solvers yet")
	else
		local SelfGrav = require 'hydro.op.selfgrav'
		self.gravOp = SelfGrav{solver = self.solver}
		self.solver.ops:insert(self.gravOp)

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

function Euler:createInitState()
	Euler.super.createInitState(self)
	local double = false --solver.app.real == 'double'
	self:addGuiVars{	
		{name='heatCapacityRatio', value=7/5},				-- unitless
		{name='rhoMin', value=double and 1e-15 or 1e-7, units='kg/m^3'},
		{name='PMin', value=double and 1e-15 or 1e-7, units='kg/(m*s^2)'},
	}
end

function Euler:initCodeModules()
	Euler.super.initCodeModules(self)

	-- added by request only, so I don't have to compile the real3x3 code
	-- not used at the moment
	self.solver.modules:add{
		name = 'calcCellMinMaxEigenvalues',
		depends = {
			'real3x3',
			'eqn.prim-cons',
		},
		code = self:template[[
range_t calcCellMinMaxEigenvalues(
	constant <?=solver.solver_t?>* solver,
	const global <?=eqn.cons_t?>* U,
	real3 x,
	real3x3 nL,
	real3x3 nU,
	real nLen
) {
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real v_n = real3_dot(W.v, nL.x);
	real Cs = calc_Cs(solver, &W);
	real Cs_nLen = Cs * nLen;
	return (range_t){
		.min = v_n - Cs_nLen, 
		.max = v_n + Cs_nLen,
	};
}
]],
	}

	self.solver.modules:add{
		name = 'eigen_forCell',
		depends = {
			'normal_t',	-- normal_t
			'coord_lower',
			'eqn.cons_t',
			'eqn.prim_t',
			'eqn.eigen_t',
			'eqn.prim-cons',
			'eqn.solvercode',	-- calc_hTotal
		},
		code = self:template[[
// used by PLM
<?=eqn.eigen_t?> eigen_forCell(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normal_t n
) {
	<?=eqn.prim_t?> W = primFromCons(solver, U, x);
	real3 vL = coord_lower(W.v, x);
	real vSq = real3_dot(W.v, vL);
	real v_n = normal_vecDotN1(n, W.v);
	real eKin = .5 * vSq;
	real hTotal = calc_hTotal(W.rho, W.P, U.ETotal);
	real CsSq = (solver->heatCapacityRatio - 1.) * (hTotal - eKin);
	real Cs = sqrt(CsSq);
	return (<?=eqn.eigen_t?>){
		.rho = W.rho,
		.v = W.v,
		.vSq = vSq,
		.vL = vL,
		.hTotal = hTotal,
		.Cs = Cs,
	};
}
]],
	}
end

function Euler:initCodeModule_fluxFromCons()
	self.solver.modules:add{
		name = 'fluxFromCons',
		depends = {
			'solver.solver_t',
			'eqn.prim-cons',
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
	real v_n = normal_vecDotN1(n, W.v);
	real HTotal = U.ETotal + W.P;
	
	return (<?=eqn.cons_t?>){
		.rho = U.rho * v_n,
		.m = real3_add(
			real3_real_mul(U.m, v_n),
			_real3(
				normal_u1x(n) * W.P,
				normal_u1y(n) * W.P,
				normal_u1z(n) * W.P
			)
		),
		.ETotal = HTotal * v_n,
		.ePot = 0,
	};
}
]],
	}
end

function Euler:initCodeModuleCommon()
	self.solver.modules:add{
		name = 'eqn.common',
		depends = {
			'eqn.cons_t',
			'eqn.prim_t',
			'eqn.waves_t',
			'eqn.eigen_t',
			'coordLenSq',
		},
		code = self:template[[
real calc_H(constant <?=solver.solver_t?>* solver, real P) { return P * (solver->heatCapacityRatio / (solver->heatCapacityRatio - 1.)); }
real calc_h(constant <?=solver.solver_t?>* solver, real rho, real P) { return calc_H(solver, P) / rho; }
real calc_HTotal(real P, real ETotal) { return P + ETotal; }
real calc_hTotal(real rho, real P, real ETotal) { return calc_HTotal(P, ETotal) / rho; }
real calc_eKin(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.v, x); }
real calc_EKin(<?=eqn.prim_t?> W, real3 x) { return W.rho * calc_eKin(W, x); }
real calc_EInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return W.P / (solver->heatCapacityRatio - 1.); }
real calc_eInt(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W) { return calc_EInt(solver, W) / W.rho; }
real calc_EKin_fromCons(<?=eqn.cons_t?> U, real3 x) { return .5 * coordLenSq(U.m, x) / U.rho; }
real calc_ETotal(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) {
	return calc_EKin(W, x) + calc_EInt(solver, W);
}

real calc_Cs(constant <?=solver.solver_t?>* solver, const <?=eqn.prim_t?>* W) {
	return sqrt(solver->heatCapacityRatio * W->P / W->rho);
}

real calc_P(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U, real3 x) {
	real EKin = calc_EKin_fromCons(U, x);
	real EInt = U.ETotal - EKin;
	return (solver->heatCapacityRatio - 1.) * EInt;
}

]],
	}
end

-- override the default
function Euler:initCodeModulePrimCons()
	self.solver.modules:add{
		name = 'eqn.prim-cons',
		depends = {
			'real3',
			'solver.solver_t',
			'eqn.prim_t',
			'eqn.cons_t',
			'eqn.common',	-- all the calc_* stuff
		},
		code = self:template[[
<?=eqn.prim_t?> primFromCons(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U, real3 x) {
	return (<?=eqn.prim_t?>){
		.rho = U.rho,
		.v = real3_real_mul(U.m, 1./U.rho),
		.P = calc_P(solver, U, x),
		.ePot = U.ePot,
	};
}

<?=eqn.cons_t?> consFromPrim(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) {
	return (<?=eqn.cons_t?>){
		.rho = W.rho,
		.m = real3_real_mul(W.v, W.rho),
		.ETotal = calc_ETotal(solver, W, x),
		.ePot = W.ePot,
	};
}
]],
	}

	-- only used by PLM
	self.solver.modules:add{
		name = 'eqn.dU-dW',
		depends = {
			'solver.solver_t',
			'eqn.prim_t',
			'eqn.cons_t',
			'coord_lower',
		},
		code = self:template[[
<?=eqn.cons_t?> apply_dU_dW(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA, 
	<?=eqn.prim_t?> W, 
	real3 x
) {
	real3 WA_vL = coord_lower(WA.v, x);
	return (<?=eqn.cons_t?>){
		.rho = W.rho,
		.m = real3_add(
			real3_real_mul(WA.v, W.rho), 
			real3_real_mul(W.v, WA.rho)),
		.ETotal = W.rho * .5 * real3_dot(WA.v, WA_vL) 
			+ WA.rho * real3_dot(W.v, WA_vL)
			+ W.P / (solver->heatCapacityRatio - 1.),
		.ePot = W.ePot,
	};
}

<?=eqn.prim_t?> apply_dW_dU(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA,
	<?=eqn.cons_t?> U,
	real3 x
) {
	real3 WA_vL = coord_lower(WA.v, x);
	return (<?=eqn.prim_t?>){
		.rho = U.rho,
		.v = real3_sub(
			real3_real_mul(U.m, 1. / WA.rho),
			real3_real_mul(WA.v, U.rho / WA.rho)),
		.P = (solver->heatCapacityRatio - 1.) * (
			.5 * real3_dot(WA.v, WA_vL) * U.rho 
			- real3_dot(U.m, WA_vL)
			+ U.ETotal),
		.ePot = U.ePot,
	};
}
]],
	}
end


function Euler:getModuleDependsApplyInitCond()
	return table(Euler.super.getModuleDependsApplyInitCond(self))
	:append{'cartesianToCoord'}
end
Euler.initCondCode = [[
kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real3 mids = real3_real_mul(real3_add(solver->initCondMins, solver->initCondMaxs), .5);
	bool lhs = true<?
for i=1,solver.dim do
	local xi = xNames[i]
?> && x.<?=xi?> < mids.<?=xi?><?
end
?>;

	// these are all standard for all init/euler.lua initial conditions
	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real3 D = real3_zero;
	real3 B = real3_zero;
	real ePot = 0;

	<?=code?>

	<?=eqn.prim_t?> W = {
		.rho = rho,
		.v = cartesianToCoord(v, x),
		.P = P,
		.ePot = ePot,
	};

	*U = consFromPrim(solver, W, x);
}
]]

Euler.solverCodeFile = 'hydro/eqn/euler.cl'

Euler.displayVarCodeUsesPrims = true

-- [=[
Euler.predefinedDisplayVars = {
-- [[	
	'U rho',
	'U v',
	'U P',
--]]
--[[
	-- now that I've switched to components, I can't display these
	-- TODO ...unless I allow for multiple displays of the same displayVar...
	-- TODO TODO after building kernels, duplicate all display vars for each component type
	-- then things will look just like they did before, except with less kernels
	'U v x',
	'U v y',
	'U v z',
--]]
--[[	
	'U ePot',
	'U EPot',
	'U gravity',
--]]
	'U div v',
	'U curl v',
}
--]=]

function Euler:getDisplayVars()
	local vars = Euler.super.getDisplayVars(self)
	vars:append{
		-- TODO should the default display generated of variables be in solver units or SI units?
		-- should SI unit displays be auto generated as well?
		{name='v', code='value.vreal3 = W.v;', type='real3', units='m/s'},
		{name='P', code='value.vreal = W.P;', units='kg/(m*s^2)'},
		{name='eInt', code='value.vreal = calc_eInt(solver, W);', units='m^2/s^2'},
		{name='eKin', code='value.vreal = calc_eKin(W, x);', units='m^2/s^2'},
		{name='eTotal', code='value.vreal = U->ETotal / W.rho;', units='m^2/s^2'},
		{name='EInt', code='value.vreal = calc_EInt(solver, W);', units='kg/(m*s^2)'},
		{name='EKin', code='value.vreal = calc_EKin(W, x);', units='kg/(m*s^2)'},
		{name='EPot', code='value.vreal = U->rho * U->ePot;', units='kg/(m*s^2)'},
		{name='S', code='value.vreal = W.P / pow(W.rho, (real)solver->heatCapacityRatio);'},
		{name='H', code='value.vreal = calc_H(solver, W.P);', units='kg/(m*s^2)'},
		{name='h', code='value.vreal = calc_h(solver, W.rho, W.P);', units='m^2/s^2'},
		{name='HTotal', code='value.vreal = calc_HTotal(W.P, U->ETotal);', units='kg/(m*s^2)'},
		{name='hTotal', code='value.vreal = calc_hTotal(W.rho, W.P, U->ETotal);', units='m^2/s^2'},
		{name='speed of sound', code='value.vreal = calc_Cs(solver, &W);', units='m/s'},
		{name='Mach number', code='value.vreal = coordLen(W.v, x) / calc_Cs(solver, &W);'},
		{name='temperature', code=self:template[[
<? local clnumber = require 'cl.obj.number' ?>
<? local materials = require 'hydro.materials' ?>
#define C_v				<?=('%.50f'):format(materials.Air.C_v)?>
	value.vreal = calc_eInt(solver, W) / C_v;
]], units='K'},
	}:append(self.gravOp and
		{{name='gravity', code=self:template[[
	if (!OOB(1,1)) {
		value.vreal3 = calcGravityAccel<?=eqn.gravOp.name?>(solver, U);
	}
]], type='real3', units='m/s^2'}} or nil
	)

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

	return vars
end

Euler.eigenVars = table{
	-- Roe-averaged vars
	{name='rho', type='real', units='kg/m^3'},
	{name='v', type='real3', units='m/s'},
	{name='hTotal', type='real', units='m^2/s^2'},
	-- derived vars
	{name='Cs', type='real', units='m/s'},
	{name='vSq', type='real', units='m^2/s^2'},
	{name='vL', type='real3', units='m/s'},
}

function Euler:eigenWaveCodePrefix(n, eig, x)
	return self:template([[
	real Cs_nLen = <?=eig?>.Cs * normal_len(<?=n?>);
	real v_n = normal_vecDotN1(<?=n?>, <?=eig?>.v);
]],	{
		eig = '('..eig..')',
		x = x,
		n = n,
	})
end

-- W is an extra param specific to Euler's calcDT in this case
-- but then I just explicitly wrote out the calcDT, so the extra parameters just aren't used anymore.
function Euler:consWaveCodePrefix(n, U, x)
	return self:template([[
	<?=eqn.prim_t?> W = primFromCons(solver, <?=U?>, <?=x?>);
	real Cs_nLen = calc_Cs(solver, &W) * normal_len(<?=n?>);
	real v_n = normal_vecDotN1(<?=n?>, W.v);
]], {
		U = '('..U..')',
		n = n,
		x = x,
	})
end

function Euler:consWaveCode(n, U, x, waveIndex)
	if waveIndex == 0 then
		return '(v_n - Cs_nLen)'
	elseif waveIndex >= 1 and waveIndex <= 3 then
		return 'v_n'
	elseif waveIndex == 4 then
		return '(v_n + Cs_nLen)'
	end
	error'got a bad waveIndex'
end

-- as long as U or eig isn't used, we can use this for both implementations
Euler.eigenWaveCode = Euler.consWaveCode

-- this one calcs cell prims once and uses it for all sides
-- it is put here instead of in hydro/eqn/euler.cl so euler-burgers can override it
-- TODO move the sqrt() out of the loop altogether?
function Euler:initCodeModule_calcDT()
	local solver = self.solver
	solver.modules:add{
		name = 'calcDT',
		depends = table{
			'solver.solver_t',
			'eqn.prim-cons',
			'eqn.guiVars.compileTime',
		}:append(
			require 'hydro.solver.meshsolver'.is(solver) and {
				'normal_t',
			} or nil
		),
		code = self:template[[
<? if require 'hydro.solver.gridsolver'.is(solver) then ?>

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
	real3 x = cellBuf[index].pos;

	const global <?=eqn.cons_t?>* U = UBuf + index;
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real Cs = calc_Cs(solver, &W);

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
<? 
if solver.coord.vectorComponent == 'cartesian' 
and not require 'hydro.coord.cartesian'.is(solver.coord)
then 
?>		real dx = cell_dx<?=side?>(x); 
<? else 
?>		real dx = solver->grid_dx.s<?=side?>;
<? end 
?>
		if (dx > 1e-7) {
			//use cell-centered eigenvalues
			real v_n = normal_vecDotN1(normal_forSide<?=side?>(x), W.v);
			real lambdaMin = v_n - Cs;
			real lambdaMax = v_n + Cs;
			real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
			absLambdaMax = max((real)1e-9, absLambdaMax);
			//TODO this should be based on coord + vectorComponent 
			//non-cartesian coord + cartesian component uses |u(x+dx)-u(x)|
			dt = (real)min(dt, dx / absLambdaMax);
		}
	}<? end ?>
	dtBuf[index] = dt;
}

<? else -- mesh solver ?>

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,					//[numCells]
	const global <?=eqn.cons_t?>* UBuf,	//[numCells]
	const global <?=solver.coord.cell_t?>* cellBuf,		//[numCells]
	const global <?=solver.coord.face_t?>* faces,			//[numFaces]
	const global int* cellFaceIndexes	//[numCellFaceIndexes]
) {
	int cellIndex = get_global_id(0);
	if (cellIndex >= get_global_size(0)) return;
	
	const global <?=eqn.cons_t?>* U = UBuf + cellIndex;
	const global <?=solver.coord.cell_t?>* cell = cellBuf + cellIndex;
	real3 x = cell->pos;

	real dt = INFINITY;
	for (int i = 0; i < cell->faceCount; ++i) {
		const global <?=solver.coord.face_t?>* face = faces + cellFaceIndexes[i + cell->faceOffset];
		real dx = face->area;
		if (dx > 1e-7 && face->cells.x != -1 && face->cells.y != -1) {
			//all sides? or only the most prominent side?
			//which should we pick eigenvalues from?
			//use cell-centered eigenvalues
			normal_t n = normal_forFace(face);
			<?=eqn:consWaveCodePrefix('n', '*U', 'x')?>
			real lambdaMin = <?=eqn:consMinWaveCode('n', '*U', 'x')?>;
			real lambdaMax = <?=eqn:consMaxWaveCode('n', '*U', 'x')?>;
			real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
			absLambdaMax = max((real)1e-9, absLambdaMax);
			dt = (real)min(dt, dx / absLambdaMax);
		}
	}
	dtBuf[cellIndex] = dt;
}


<? end -- mesh vs grid solver ?>
]],
	}
end

function Euler:getModuleDependsSolver() 
	return table{
		'eqn.prim-cons',
		'normal_t',
		'coord_lower',
	}:append(
		not require 'hydro.coord.cartesian'.is(self.solver.coord)
		and	{
			'coord_conn_apply13',
			'coord_conn_apply23',
			'coord_conn_trace23',
		} or nil
	)
end

return Euler

--[[
scratch work:

P V = n R T = m R_spec T
ideal gas EOS: P V_m = R T
T = T_C + 273.15' C
R = gas constant = 8.314459848 J / (mol K) 
R_spec = 287.058 J / (kg K) for dry air 
R_spec = R / M = k_B / m
M = molar mass
k_B = Boltzmann constant = 1.3806485279e-23 J / K
m = mass within volume V	<-> ρ = m / V
R_spec = k_B ρ V

caloric perfect:
P = (gamma - 1) ρ e_int
gamma = C_p / C_v = adiabatic index / ratio of specific heats
e_int = C_v T = internal energy per unit mass / specific internal energy
C_v = specific heat at constant volume
C_p = specific heat at constant pressure
P = ρ R_spec T
e_int = C_v T
P = ρ (C_p - C_v) e_int / C_v = ρ (C_p / C_v - 1) e_int

0 C = 273.15 K

using some real-world numbers ...
P = (gamma - 1) ρ e_int
101325 kg / (m s^2) = (C_p - C_v) / C_v (1.2754 kg / m^3) C_v T 
101325 kg / (m s^2) = (1006 - 717.1) J / (kg K) (1.2754 kg / m^3) T 
T = 274.99364522457 K
T = 1.8436452245715 C ... should be
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'
local xNames = require 'common'.xNames


local Euler = class(Equation)
Euler.name = 'Euler'

-- ePot is the 6th param
-- which means it's now in the derivBuf, but it is always zero
-- so TODO a new variable for deriv size vs cons_t size?
Euler.numStates = 6	

Euler.numWaves = 5
Euler.numIntStates = 5	-- don't bother integrate ePot

Euler.hasCalcDTCode = true
Euler.hasEigenCode = true
Euler.hasFluxFromConsCode = true
Euler.roeUseFluxFromCons = true

-- the only source term that the Euler equations has is the connection coefficients of the velocity vector
-- maybe later I will automatically flag what elements are vectors
-- and automatically add connection coefficients
--
-- I'm also using the source term for the viscosity ... if you choose to use FANS
Euler.useSourceTerm = true
Euler.useConstrainU = true

Euler.initStates = require 'init.euler'


-- TODO primVars doesn't autogen displayVars, and therefore units doesn't matter
Euler.primVars = {
	{name='rho', type='real', units='kg/m^3'},
	{name='v', type='real3', units='m/s', variance='u'},			-- contravariant
	{name='P', type='real', units='kg/(m*s^2)'},
	{name='ePot', type='real', units='m^2/s^2'},
}

Euler.consVars = {
	{name='rho', type='real', units='kg/m^3'},
	{name='m', type='real3', units='kg/(m^2*s)', variance='u'},	-- contravariant
	{name='ETotal', type='real', units='kg/(m*s^2)'},
	{name='ePot', type='real', units='m^2/s^2'},
}


function Euler:init(args)
	Euler.super.init(self, args)

	local SelfGrav = require 'op.selfgrav'

	if require 'solver.meshsolver'.is(self.solver) then
		print("not using selfgrav with mesh solvers yet")
	else
		self.gravOp = SelfGrav{solver = self.solver}
		self.solver.ops:insert(self.gravOp)
	end
end

function Euler:createInitState()
	Euler.super.createInitState(self)
	local double = false --solver.app.real == 'double'
	self:addGuiVars{	
		{name='heatCapacityRatio', value=7/5},				-- unitless
		{name='rhoMin', value=double and 1e-15 or 1e-7},	-- kg/m^3
		{name='PMin', value=double and 1e-15 or 1e-7},		-- kg/m^3
	}
end

function Euler:getCommonFuncCode()
	return template([[
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

]], {
		solver = self.solver,
		eqn = self,
	})
end

function Euler:getPrimConsCode()
	return template([[

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


]], {
		solver = self.solver,
		eqn = self,
	})
end

Euler.initStateCode = [[
<? local xNames = require 'common'.xNames ?>
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
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

	UBuf[index] = consFromPrim(solver, W, x);
}
]]

Euler.solverCodeFile = 'eqn/euler.cl'

Euler.displayVarCodeUsesPrims = true

-- k is 0,1,2
local function vorticity(eqn,k,result)
	local i = (k+1)%3
	local j = (i+1)%3
	return {
		name = 'vorticity '..xNames[k+1],
		code = template([[
	if (OOB(1,1)) {
		<?=result?> = 0.;
	} else {
		global const <?=eqn.cons_t?>* Uim = U - solver->stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Uip = U + solver->stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Ujm = U - solver->stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + solver->stepsize.s<?=j?>;

		//TODO incorporate metric

		real vim_j = Uim->m.s<?=j?> / Uim->rho;
		real vip_j = Uip->m.s<?=j?> / Uip->rho;
		
		real vjm_i = Ujm->m.s<?=i?> / Ujm->rho;
		real vjp_i = Ujp->m.s<?=i?> / Ujp->rho;
		
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

-- [=[
Euler.predefinedDisplayVars = {
-- [[	
	'U rho',
	'U v',
	'U P',
--]]
	'U v mag',
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
	}:append(self.gravOp and
		{{name='gravity', code=template([[
	if (!OOB(1,1)) {
		value.vreal3 = calcGravityAccel<?=eqn.gravOp.name?>(solver, U);
	}
]], {eqn=self}), type='real3', units='m/s'}} or nil
	):append{
		{name='temp', code=template([[
<? local clnumber = require 'cl.obj.number' ?>
<? local materials = require 'materials' ?>
#define C_v				<?=('%.50f'):format(materials.Air.C_v)?>
	value.vreal = calc_eInt(solver, W) / C_v;
]]), units='K'}
	}

	-- vorticity = [x ,y ,z] [v.x, v.y, v.z][
	-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
	
	if not require 'solver.meshsolver'.is(self.solver) then
		if self.solver.dim == 2 then
			vars:insert(vorticity(self,2,'value.vreal'))
		elseif self.solver.dim == 3 then
			local v = xNames:mapi(function(xi,i)
				return vorticity(self,i-1,'value.vreal3.'..xi)
			end)
			vars:insert{name='vorticityVec', code=template([[
	<? for i,vi in ipairs(v) do ?>{
		<?=vi.code?>
	}<? end ?>
]], {v=v}), type='real3', units='m/s^2'}
		end
	end

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

function Euler:eigenWaveCodePrefix(n, eig, x, W)
	return template([[
	real Cs_nLen = <?=eig?>.Cs * normalInfo_len(<?=n?>);
	real v_n = normalInfo_vecDotN1(<?=n?>, <?=eig?>.v);
]],	{
		eqn = self,
		eig = '('..eig..')',
		W = W and '('..W..')' or nil,
		x = x,
		n = n,
	})
end

-- W is an extra param specific to Euler's calcDT in this case
function Euler:consWaveCodePrefix(n, U, x, W)
	return template([[
<? if not W then 
	W = 'W'
?>
	<?=eqn.prim_t?> W = primFromCons(solver, <?=U?>, <?=x?>);
<? end ?>
	real Cs_nLen = calc_Cs(solver, &<?=W?>) * normalInfo_len(<?=n?>);
	real v_n = normalInfo_vecDotN1(<?=n?>, <?=W?>.v);
]], {
		eqn = self,
		U = '('..U..')',
		W = W and '('..W..')' or nil,
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
-- it is put here instead of in eqn/euler.cl so euler-burgers can override it
-- TODO move the sqrt() out of the loop altogether?
function Euler:getCalcDTCode()
	return template([[
<? local solver = eqn.solver ?>
<? if require 'solver.gridsolver'.is(solver) then ?>

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 x = cell_x(i);

	const global <?=eqn.cons_t?>* U = UBuf + index;
	<?=eqn.prim_t?> W = primFromCons(solver, *U, x);
	real Cs = calc_Cs(solver, &W);

	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
		//use cell-centered eigenvalues
		real v_n = normalInfo_vecDotN1(normalInfo_forSide<?=side?>(x), W.v);
		real lambdaMin = v_n - Cs;
		real lambdaMax = v_n + Cs;
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
		absLambdaMax = max((real)1e-9, absLambdaMax);
		//TODO this should be based on coord + vectorComponent 
		//non-cartesian coord + cartesian component uses |u(x+dx)-u(x)|
		real dx = solver->grid_dx.s<?=side?>;
		dt = (real)min(dt, dx / absLambdaMax);
	}<? end ?>
	dtBuf[index] = dt;
}

<? else -- mesh solver ?>

kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,					//[numCells]
	const global <?=eqn.cons_t?>* UBuf,	//[numCells]
	const global cell_t* cells,			//[numCells]
	const global face_t* faces			//[numFaces]
) {
	int cellIndex = get_global_id(0);
	if (cellIndex >= get_global_size(0)) return;
	
	const global <?=eqn.cons_t?>* U = UBuf + cellIndex;
	const global cell_t* cell = cells + cellIndex;

	real dt = INFINITY;
	for (int i = 0; i < cell->numSides; ++i) {
		const global face_t* face = faces + cell->faces[i];
		//all sides? or only the most prominent side?
		//which should we pick eigenvalues from?
		//use cell-centered eigenvalues
		real3 x = cell->x;
		normalInfo_t n = normalInfo_forCell(cell);
		<?=eqn:consWaveCodePrefix('n', '*U', 'x')?>
		real lambdaMin = <?=eqn:consMinWaveCode('n', '*U', 'x')?>;
		real lambdaMax = <?=eqn:consMaxWaveCode('n', '*U', 'x')?>;
		real absLambdaMax = max(fabs(lambdaMin), fabs(lambdaMax));
		absLambdaMax = max((real)1e-9, absLambdaMax);
		real dx = cell->maxDist;
		dt = (real)min(dt, dx / absLambdaMax);
	}
	dtBuf[cellIndex] = dt;
}


<? end -- mesh vs grid solver ?>
]], {
		eqn = self,
	})
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
m = mass within volume V	<-> rho = m / V
R_spec = k_B rho V

caloric perfect:
P = (gamma - 1) rho e_int
gamma = C_p / C_v = adiabatic index / ratio of specific heats
e_int = C_v T = internal energy per unit mass / specific internal energy
C_v = specific heat at constant volume
C_p = specific heat at constant pressure
P = rho R_spec T
e_int = C_v T
P = rho (C_p - C_v) e_int / C_v = rho (C_p / C_v - 1) e_int

0 C = 273.15 K

using some real-world numbers ...
P = (gamma - 1) rho e_int
101325 kg / (m s^2) = (C_p - C_v) / C_v (1.2754 kg / m^3) C_v T 
101325 kg / (m s^2) = (1006 - 717.1) J / (kg K) (1.2754 kg / m^3) T 
T = 274.99364522457 K
T = 1.8436452245715 C ... should be
--]]

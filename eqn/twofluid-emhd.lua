--[[
2014 Abgrall, Kumar "Robust Finite Volume Schemes for Two-Fluid Plasma Equations"

the solver/twofluid-emhd-behavior.lua instanciates three separate solvers:
	ion, electron, maxwell
However OpenCL is having trouble with this.
Maybe I have an oob memory write somewhere?
Maybe running code from two separate programs doesn't propery block on the GPU?
Either way, here's the three equations combined into one.
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local clnumber = require 'cl.obj.number'
local template = require 'template'

local fluids = table{'ion', 'elec'}

local TwoFluidEMHD = class(Equation)

-- set this to 'true' to use the init.euler states
-- these states are only provided in terms of a single density and pressure variable,
-- so the subsequent ion and electron densities and pressures must be derived from this.
TwoFluidEMHD.useEulerInitState = true

TwoFluidEMHD.name = 'TwoFluidEMHD'
TwoFluidEMHD.numWaves = 18
TwoFluidEMHD.numIntStates = 18

TwoFluidEMHD.consVars = table{
	--integration variables		
	{ion_rho = 'real'},
	{ion_m = 'real3'},
	{ion_ETotal = 'real'},
	
	{elec_rho = 'real'},
	{elec_m = 'real3'},
	{elec_ETotal = 'real'},

	{E = 'real3'},
	{B = 'real3'},
	{phi = 'real'},	-- E potential
	{psi = 'real'},	-- B potential

	--extra	
	{ion_ePot = 'real'},
	{elec_ePot = 'real'},
}

TwoFluidEMHD.primVars = table{
	--integration variables		
	{ion_rho = 'real'},
	{ion_v = 'real3'},
	{ion_P = 'real'},
	
	{elec_rho = 'real'},
	{elec_v = 'real3'},
	{elec_P = 'real'},

	{B = 'real3'},
	{E = 'real3'},
	{phi = 'real'},
	{psi = 'real'},
	
	--extra	
	{ion_ePot = 'real'},
	{elec_ePot = 'real'},
}

TwoFluidEMHD.mirrorVars = {
	{'ion_m.x', 'elec_m.x', 'E.x', 'B.x'}, 
	{'ion_m.y', 'elec_m.y', 'E.y', 'B.y'}, 
	{'ion_m.z', 'elec_m.z', 'E.z', 'B.z'},
}

TwoFluidEMHD.hasEigenCode = true
TwoFluidEMHD.useSourceTerm = true
TwoFluidEMHD.useConstrainU = true
TwoFluidEMHD.roeUseFluxFromCons = true

if TwoFluidEMHD.useEulerInitState then
	TwoFluidEMHD.initStates = require 'init.euler'
else
	TwoFluidEMHD.initStates = require 'init.twofluid-emhd'
end

-- TODO this has the symptoms of all the other CL kernels that intel's compiler bugged out on
-- i.e. takes a minute to run a kernel on the GPU
function TwoFluidEMHD:init(solver)
	TwoFluidEMHD.super.init(self, solver)

--	local NoDiv = require 'solver.nodiv'
--	solver.ops:insert(NoDiv{solver=solver})	-- nodiv on maxwell ... or just use potentials 

io.stderr:write'you need to give different selfgravs different names for twofluid selfgrav to work\n' 
	--[[ TODO give each selfgrav a unique function name
	local SelfGrav = require 'solver.selfgrav'
	solver.ops:insert(SelfGrav{
		solver = solver,
		densityField = 'ion_rho',
		potentialField = 'ion_ePot',
	})	-- selfgrav on ion 
	solver.ops:insert(SelfGrav{
		solver = solver,
		densityField = 'elec_rho',
		potentialField = 'elec_ePot',
	})	-- selfgrav on electron
	--]]
end

function TwoFluidEMHD:createInitState()
	TwoFluidEMHD.super.createInitState(self)
	
	self:addGuiVars(table{
		--never given, only stated as "speeds for the Maxwell equation"
		-- of course, they're associated with the potentials, so... they could be arbitrary
		{name='divPsiWavespeed', value=1},
		{name='divPhiWavespeed', value=1},
		
		-- gamma = heat capacity ratio
		{name='heatCapacityRatio', value=5/3},

		-- x_0 = reference length
		{name='referenceLength', value=1},

		-- B_0 = reference magnetic field
		{name='referenceMagneticField', value=1},

		-- v_i^T = ion reference thermal velocity
		--{name='ionReferenceThermalVelocity', value=1},
		-- normalized Larmor radius:
		-- lHat = l_r / x_0 = m_i v_i^T / (q_i B_0 x_0)	
		-- therefore:
		-- v_i^T = l_r B_0 q_i / m_i
		-- just use this equation:
		-- v_i^T = l_r r_i B_0

		-- m_i = ion mass
		{name='ionMass', value=1},
	
		-- c = speed of light
		{name='speedOfLight', value=1},
		-- normalized speed of light:
		-- cHat = c / v_i^T
		
		-- l_r = larmor radius
		{name='ionLarmorRadius', value=1},
		
		-- m = m_i / m_e
		-- https://en.wikipedia.org/wiki/Proton-to-electron_mass_ratio
		-- m_i / m_e = 1836.15267389
		{name='ionElectronMassRatio', value=100},
		
		-- r_i = q_i / m_i
		{name='ionChargeMassRatio', value=1},
		
		-- r_e = q_e / m_e
		-- https://en.wikipedia.org/wiki/Mass-to-charge_ratio
		-- q_e / m_e = -1.758820024e+11 C/kg
		{name='elecChargeMassRatio', value=.05},
	
		-- lambda_d = ion Debye length
		{name='ionDebyeLength', value=1},
		-- normalized ion Debye length: 
		-- lambdaHat_d = lambda_d / l_r
	}:append(fluids:map(function(fluid)
		return table{
			{name='min_'..fluid..'_rho', value=1e-4},
			{name='min_'..fluid..'_P', value=1e-4},
		}
	end):unpack()))
end

function TwoFluidEMHD:getCommonFuncCode()
	return template([[
real ESq(<?=eqn.cons_t?> U, real3 x) { return real3_lenSq(U.E); }
real BSq(<?=eqn.cons_t?> U, real3 x) { return real3_lenSq(U.B); }

inline real calc_H(real P) { return P * (heatCapacityRatio / (heatCapacityRatio - 1.)); }
inline real calc_h(real rho, real P) { return calc_H(P) / rho; }
inline real calc_hTotal(real rho, real P, real ETotal) { return (P + ETotal) / rho; }
inline real calc_HTotal(real P, real ETotal) { return P + ETotal; }

<? for _,fluid in ipairs(fluids) do ?>
inline real calc_<?=fluid?>_eKin(<?=eqn.prim_t?> W, real3 x) { return .5 * coordLenSq(W.<?=fluid?>_v, x); }
inline real calc_<?=fluid?>_EKin(<?=eqn.prim_t?> W, real3 x) { return W.<?=fluid?>_rho * calc_<?=fluid?>_eKin(W, x); }
inline real calc_<?=fluid?>_EInt(<?=eqn.prim_t?> W) { return W.<?=fluid?>_P / (heatCapacityRatio - 1.); }
inline real calc_<?=fluid?>_eInt(<?=eqn.prim_t?> W) { return calc_<?=fluid?>_EInt(W) / W.<?=fluid?>_rho; }
inline real calc_<?=fluid?>_EKin_fromCons(<?=eqn.cons_t?> U, real3 x) { return .5 * coordLenSq(U.<?=fluid?>_m, x) / U.<?=fluid?>_rho; }
inline real calc_<?=fluid?>_ETotal(<?=eqn.prim_t?> W, real3 x) {
	real EPot = W.<?=fluid?>_rho * W.<?=fluid?>_ePot;
	return calc_<?=fluid?>_EKin(W, x) + calc_<?=fluid?>_EInt(W) + EPot;
}
inline real calc_<?=fluid?>_Cs(const <?=eqn.prim_t?>* W) {
	return sqrt(heatCapacityRatio * W-><?=fluid?>_P / W-><?=fluid?>_rho);
}
<? end ?>
]], {
		eqn = self,
		fluids = fluids,
	})
end

function TwoFluidEMHD:getPrimConsCode()
	return template([[
inline <?=eqn.prim_t?> primFromCons(<?=eqn.cons_t?> U, real3 x) {
	<? for _,fluid in ipairs(fluids) do ?>
	real <?=fluid?>_EPot = U.<?=fluid?>_rho * U.<?=fluid?>_ePot;
	real <?=fluid?>_EKin = calc_<?=fluid?>_EKin_fromCons(U, x);
	real <?=fluid?>_EInt = U.<?=fluid?>_ETotal - <?=fluid?>_EKin - <?=fluid?>_EPot;
	<? end ?>
	return (<?=eqn.prim_t?>){
		<? for _,fluid in ipairs(fluids) do ?>
		.<?=fluid?>_rho = U.<?=fluid?>_rho,
		.<?=fluid?>_v = real3_scale(U.<?=fluid?>_m, 1./U.<?=fluid?>_rho),
		.<?=fluid?>_P = (heatCapacityRatio - 1.) * <?=fluid?>_EInt,
		.<?=fluid?>_ePot = U.<?=fluid?>_ePot,
		<? end ?>
		.E = U.E,
		.B = U.B,
		.psi = U.psi,
		.phi = U.phi,
	};
}

inline <?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> W, real3 x) {
	return (<?=eqn.cons_t?>){
		<? for _,fluid in ipairs(fluids) do ?>
		.<?=fluid?>_rho = W.<?=fluid?>_rho,
		.<?=fluid?>_m = real3_scale(W.<?=fluid?>_v, W.<?=fluid?>_rho),
		.<?=fluid?>_ETotal = calc_<?=fluid?>_ETotal(W, x),
		.<?=fluid?>_ePot = W.<?=fluid?>_ePot,
		<? end ?>
		.E = W.E,
		.B = W.B,
		.psi = W.psi,
		.phi = W.phi,
	};
}
]], {
		eqn = self,
		fluids = fluids,
	})
end

-- overridden because it adds some extra parameters to the template args
-- should I either make a function for the template arg params
-- or maybe I shouldn't have super-class'd the initState code to begin with ...
function TwoFluidEMHD:getInitStateCode()
	local code = self.initState:initState(self.solver)	
	return template([[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	bool lhs = x.x < mids.x
#if dim > 1
		&& x.y < mids.y
#endif
#if dim > 2
		&& x.z < mids.z
#endif
	;
<? 
if eqn.useEulerInitState then 
?>
	real rho = 0.;
	real3 v = _real3(0,0,0);
	real P = 0;
	real ePot = 0;
<?
else
	 for _,fluid in ipairs(fluids) do
?>	real <?=fluid?>_rho = 0;
	real3 <?=fluid?>_v = _real3(0,0,0);
	real <?=fluid?>_P = 0;
	real <?=fluid?>_ePot = 0;
<? 
	end 
end
?>	real3 E = _real3(0,0,0);
	real3 B = _real3(0,0,0);
	real conductivity = 1.;
	real permittivity = 1. / (4. * M_PI);
	real permeability = 4. * M_PI;

]]..code..[[

	<?=eqn.prim_t?> W = {
<? 
if eqn.useEulerInitState then 
?>
		.ion_rho = rho,
		.elec_rho = rho / ionElectronMassRatio, 

		// "the electron pressure is taken to be elec_P = 5 ion_rho"
		// is that arbitrary?
		.elec_P = 5. * rho,
		
		// "the ion pressure is 1/100th the electron pressure"
		// is that from the mass ratio of ion/electron?
		.ion_P = P / ionElectronMassRatio, 

		.ion_v = v,
		.elec_v = v,

		.E = E,
		.B = B,
		.phi = 0,
		.psi = 0,
		
		.ion_ePot = 0,
		.elec_ePot = 0,

<?	
else	-- expect the initState to explicitly provide the ion_ and elec_ Euler fluid variables
	for _,fluid in ipairs(fluids) do ?>
		.<?=fluid?>_rho = <?=fluid?>_rho,
		.<?=fluid?>_v = cartesianToCoord(<?=fluid?>_v, x),
		.<?=fluid?>_P = <?=fluid?>_P,
		.<?=fluid?>_ePot = <?=fluid?>_ePot,
<?
	end
end
?>
		.E = E,
		.B = B,
		.psi = 0,
	};
	UBuf[index] = consFromPrim(W, x);
}
]], {
		eqn = self,
		fluids = fluids,
		clnumber = clnumber,
	})
end

function TwoFluidEMHD:getSolverCode()
	return template(file['eqn/twofluid-emhd.cl'], {
		eqn = self, 
		solver = self.solver,
		fluids = fluids,
		clnumber = clnumber,
	})
end

function TwoFluidEMHD:getDisplayVarCodePrefix()
	return template([[
	global const <?=eqn.cons_t?>* U = buf + index;
	<?=eqn.prim_t?> W = primFromCons(*U, x);
]], {
		eqn = self,
	})
end

function TwoFluidEMHD:getDisplayVars()
	local vars = TwoFluidEMHD.super.getDisplayVars(self)

	for _,fluid in ipairs(fluids) do
	
		-- k is 0,1,2
		local function vorticity(k)
			local xs = {'x','y','z'}
			local i = (k+1)%3
			local j = (i+1)%3
			return {[fluid..' vorticity '..xs[k+1]] = template([[
	if (OOB(1,1)) {
		*value = 0.;
	} else {
		global const <?=eqn.cons_t?>* Uim = U - stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Uip = U + stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Ujm = U - stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + stepsize.s<?=j?>;

		//TODO incorporate metric

		real3 vim_j = Uim-><?=fluid?>_m.s<?=j?> / Uim-><?=fluid?>_rho;
		real3 vip_j = Uip-><?=fluid?>_m.s<?=j?> / Uip-><?=fluid?>_rho;
		real3 vjm_i = Ujm-><?=fluid?>_m.s<?=i?> / Ujm-><?=fluid?>_rho;
		real3 vjp_i = Ujp-><?=fluid?>_m.s<?=i?> / Ujp-><?=fluid?>_rho;
		
		*value = (vjp_i - vjm_i) / (2. * grid_dx<?=i?>)
				- (vip_j - vim_j) / (2. * grid_dx<?=j?>);
	}
]], 		{
				i = i,
				j = j,
				eqn = self,
				fluid = fluid,
			})}
		end
		
		vars:append{
			{[fluid..' v'] = '*valuevec = W.'..fluid..'_v;', type='real3'},
			{[fluid..' P'] = '*value = W.'..fluid..'_P;'},
			{[fluid..' eInt'] = '*value = calc_'..fluid..'_eInt(W);'},
			{[fluid..' eKin'] = '*value = calc_'..fluid..'_eKin(W, x);'},
			{[fluid..' ePot'] = '*value = U->'..fluid..'_ePot;'},
			{[fluid..' EInt'] = '*value = calc_'..fluid..'_EInt(W);'},
			{[fluid..' EKin'] = '*value = calc_'..fluid..'_EKin(W, x);'},
			{[fluid..' EPot'] = '*value = U->'..fluid..'_rho * U->'..fluid..'_ePot;'},
			{[fluid..' ETotal'] = '*value = U->'..fluid..'_ETotal;'},
			{[fluid..' S'] = '*value = W.'..fluid..'_P / pow(W.'..fluid..'_rho, (real)heatCapacityRatio);'},
			{[fluid..' H'] = '*value = calc_H(W.'..fluid..'_P);'},
			{[fluid..' h'] = '*value = calc_h(W.'..fluid..'_rho, W.'..fluid..'_P);'},
			{[fluid..' HTotal'] = '*value = calc_HTotal(W.'..fluid..'_P, U->'..fluid..'_ETotal);'},
			{[fluid..' hTotal'] = '*value = calc_hTotal(W.'..fluid..'_rho, W.'..fluid..'_P, U->'..fluid..'_ETotal);'},
			{[fluid..'Speed of Sound'] = '*value = calc_'..fluid..'_Cs(&W);'},
			{[fluid..'Mach number'] = '*value = coordLen(W.'..fluid..'_v, x) / calc_'..fluid..'_Cs(&W);'},
		}:append( ({
		-- vorticity = [,x ,y ,z] [v.x, v.y, v.z][
		-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
				[1] = {},
				[2] = {vorticity(2)},
				[3] = range(0,2):map(vorticity),

		})[self.solver.dim] )
	end

	vars:append{
		{['EM energy'] = [[
	//*value = .5 * (coordLen(U->E) + coordLen(U->B));
	*value = .5 * (real3_len(U->E) + real3_len(U->B));
]]},
	}:append(table{'E','B'}:map(function(var,i)
		local field = assert( ({E='E', B='B'})[var] )
		return {['div '..var] = template([[
	*value = .5 * (0.
<?
for j=0,solver.dim-1 do
?>		+ (U[stepsize.s<?=j?>].<?=field?>.s<?=j?> 
			- U[-stepsize.s<?=j?>].<?=field?>.s<?=j?>
		) / grid_dx<?=j?>
<?
end 
?>	);
]], {solver=self.solver, field=field})}
	end))

	return vars
end

local eigenVars = table()
for _,fluid in ipairs(fluids) do
	eigenVars:append{
		-- Roe-averaged vars
		{[fluid..'_rho'] = 'real'},
		{[fluid..'_v'] = 'real3'},
		{[fluid..'_hTotal'] = 'real'},

		-- derived vars
		{[fluid..'_vSq'] = 'real'},
		{[fluid..'_Cs'] = 'real'},
	}
end

TwoFluidEMHD.eigenVars = eigenVars

function TwoFluidEMHD:eigenWaveCodePrefix(side, eig, x)
	return template([[
<? for i,fluid in ipairs(fluids) do ?>
	real <?=fluid?>_Cs_sqrt_gU = <?=eig?>-><?=fluid?>_Cs * coord_sqrt_gU<?=side..side?>(x);
	real <?=fluid?>_v_n = <?=eig?>-><?=fluid?>_v.s[<?=side?>];
<? end ?>
]], {
		eig = '('..eig..')',
		fluids = fluids,
		side = side,
	})
end

function TwoFluidEMHD:eigenWaveCode(side, eig, x, waveIndex)
	for i,fluid in ipairs(fluids) do
		if waveIndex == 0 + 5 * (i-1) then
			return template('<?=fluid?>_v_n - <?=fluid?>_Cs_sqrt_gU', {fluid=fluid})
		elseif waveIndex >= 1 + 5 * (i-1)
		and waveIndex <= 3 + 5 * (i-1)
		then
			return template('<?=fluid?>_v_n', {fluid=fluid})
		elseif waveIndex == 4 + 5 * (i-1) then
			return template('<?=fluid?>_v_n + <?=fluid?>_Cs_sqrt_gU', {fluid=fluid})
		end
	end
	if waveIndex >= 5*#fluids and waveIndex < 5*#fluids+8 then
		-- 2014 Abgrall, Kumar eqn 1.9 says the eigenvalues are c, while the flux contains cHat ...
		return ({
			'-normalizedSpeedOfLight * divPhiWavespeed',
			'-normalizedSpeedOfLight * divPsiWavespeed',
			'-normalizedSpeedOfLight',
			'-normalizedSpeedOfLight',
			'normalizedSpeedOfLight',
			'normalizedSpeedOfLight',
			'normalizedSpeedOfLight * divPsiWavespeed',
			'normalizedSpeedOfLight * divPhiWavespeed',	
		})[waveIndex - 5*#fluids + 1]
	end
	error('got a bad waveIndex: '..waveIndex)
end

function TwoFluidEMHD:getFluxFromConsCode()
	return template([[
<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x
) {
	<?=eqn.prim_t?> W = primFromCons(U, x);
	<?=eqn.cons_t?> F;

<? 
for _,fluid in ipairs(fluids) do
?>	real <?=fluid?>_vj = W.<?=fluid?>_v.s<?=side?>;
	real <?=fluid?>_HTotal = U.<?=fluid?>_ETotal + W.<?=fluid?>_P;
	
	F.<?=fluid?>_rho = U.<?=fluid?>_m.s<?=side?>;
	F.<?=fluid?>_m = real3_scale(U.<?=fluid?>_m, <?=fluid?>_vj);
<? 	for i=0,2 do
?>	F.<?=fluid?>_m.s<?=i?> += coord_gU<?=i?><?=side?>(x) * W.<?=fluid?>_P;
<? 	end
?>	F.<?=fluid?>_ETotal = <?=fluid?>_HTotal * <?=fluid?>_vj;
	F.<?=fluid?>_ePot = 0.;
<? 
end
?>

	//taken from glm-maxwell instead of the 2014 Abgrall, Kumar
	real3 B = U.B;
	real3 E = U.E;
	<? if side == 0 then ?>
	F.E = _real3(normalizedSpeedOfLight * divPhiWavespeed * U.phi, normalizedSpeedOfLight * B.z, -normalizedSpeedOfLight * B.y);
	F.B = _real3(divPsiWavespeed * U.psi, -E.z, E.y);
	<? elseif side == 1 then ?>
	F.E = _real3(-normalizedSpeedOfLight * B.z, normalizedSpeedOfLight * divPhiWavespeed * U.phi, normalizedSpeedOfLight * B.x);
	F.B = _real3(E.z, divPsiWavespeed * U.psi, -E.x);
	<? elseif side == 2 then ?>
	F.E = _real3(normalizedSpeedOfLight * B.y, -normalizedSpeedOfLight * B.x, normalizedSpeedOfLight * divPhiWavespeed * U.phi);
	F.B = _real3(-E.y, E.x, divPsiWavespeed * U.psi);
	<? end ?>
	F.phi = E.s<?=side?> * divPhiWavespeed;
	F.psi = B.s<?=side?> * divPsiWavespeed * normalizedSpeedOfLight;

	F.ion_ePot = 0;
	F.elec_ePot = 0;

	return F;
}
<? end ?>
]], {
		eqn = self,
		solver = self.solver,
		fluids = fluids,
	})
end

return TwoFluidEMHD

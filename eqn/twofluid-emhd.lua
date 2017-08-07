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

local fluids = {'ion', 'elec'}

local TwoFluidEMHD = class(Equation)

-- set this to 'true' to use the init.euler states
-- these states are only provided in terms of a single density and pressure variable,
-- so the subsequent ion and electron densities and pressures must be derived from this.
TwoFluidEMHD.useEulerInitState = true

TwoFluidEMHD.name = 'TwoFluidEMHD'
TwoFluidEMHD.numStates = 22
TwoFluidEMHD.numWaves = 16
TwoFluidEMHD.numIntStates = 16
TwoFluidEMHD.mirrorVars = {
	{'ion_m.x', 'elec_m.x', 'epsE.x', 'B.x'}, 
	{'ion_m.y', 'elec_m.y', 'epsE.y', 'B.y'}, 
	{'ion_m.z', 'elec_m.z', 'epsE.z', 'B.z'},
}

-- v_i^T = ion reference thermal velocity
-- B_0 = reference magnetic field 
-- x_0 = reference length
-- m_i = ion mass = 1
TwoFluidEMHD.ionLarmorRadius = .05			-- lHat = l_r / x_0 = m_i v_i^T / (q_i B_0 x_0)
TwoFluidEMHD.ionElectronMassRatio = 5		-- m = m_i / m_e
TwoFluidEMHD.ionChargeMassRatio = 1			-- r_i = q_i / m_i
TwoFluidEMHD.elecChargeMassRatio = .05		-- r_e = q_e / m_e
TwoFluidEMHD.speedOfLight = 1				-- cHat = c / v_i^T
TwoFluidEMHD.ionDebyeLength = 1				-- lambdaHat_d = lambda_d / l_r

TwoFluidEMHD.hasEigenCode = true
TwoFluidEMHD.useSourceTerm = true
TwoFluidEMHD.hasFluxFromCons = true

if TwoFluidEMHD.useEulerInitState then
	TwoFluidEMHD.initStates = require 'init.euler'
else
	TwoFluidEMHD.initStates = require 'init.twofluid-emhd'
end

function TwoFluidEMHD:init(...)
	self.guiVars = {
		require 'guivar.float'{name='heatCapacityRatio', value=5/3}
	}
	TwoFluidEMHD.super.init(self, ...)
end


function TwoFluidEMHD:getTypeCode()
	return template([[
typedef union {
	real ptr[<?=eqn.numStates?>];
	struct {
//integration variables		
		real ion_rho;
		real3 ion_m;
		real ion_ETotal;
		
		real elec_rho;
		real3 elec_m;
		real elec_ETotal;
	
		real3 epsE;
		real3 B;
//extra	
		real ion_ePot;
		real elec_ePot;
		real BPot;
		real sigma;
		real eps;
		real mu;
	};
} <?=eqn.cons_t?>;

typedef union {
	real ptr[<?=eqn.numStates?>];
	struct {
//integration variables		
		real ion_rho;
		real3 ion_v;
		real ion_P;
		
		real elec_rho;
		real3 elec_v;
		real elec_P;
	
		real3 E;
		real3 B;
//extra	
		real ion_ePot;
		real elec_ePot;
		real BPot;
		real sigma;
		real eps;
		real mu;
	};
} <?=eqn.prim_t?>;
]], {
	eqn = self,
})
end

function TwoFluidEMHD:getCodePrefix()
	return table{
		TwoFluidEMHD.super.getCodePrefix(self),
		template([[
real ESq(<?=eqn.cons_t?> U, real3 x) { return real3_lenSq(U.epsE) / (U.eps * U.eps); }
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
		.E = real3_scale(U.epsE, 1./U.eps),
		.B = U.B,
		.BPot = U.BPot,
		.sigma = U.sigma,
		.eps = U.eps,
		.mu = U.mu,
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
		.epsE = real3_scale(W.E, W.eps),
		.B = W.B,
		.BPot = W.BPot,
		.sigma = W.sigma,
		.eps = W.eps,
		.mu = W.mu,
	};
}
]], {
	eqn = self,
	fluids = fluids,
}),
	}:concat'\n'
end

function TwoFluidEMHD:getInitStateCode()
	local initState = self.initStates[self.solver.initStateIndex]
	assert(initState, "couldn't find initState "..self.solver.initStateIndex)
	local code = initState.init(self.solver)	
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
		.elec_rho = rho * <?=clnumber(1/eqn.ionElectronMassRatio)?>, 

		// "the electron pressure is taken to be elec_P = 5 ion_rho"
		// is that arbitrary?
		.elec_P = 5. * rho,
		
		// "the ion pressure is 1/100th the electron pressure"
		// is that from the mass ratio of ion/electron?
		.ion_P = P * <?=clnumber(1/eqn.ionElectronMassRatio)?>, 

		.ion_v = v,
		.elec_v = v,

		.E = E,
		.B = B,
		.BPot = 0,
		.sigma = conductivity,
		.eps = permittivity,
		.mu = permeability,

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
		.BPot = 0,
		.sigma = conductivity,
		.eps = permittivity,
		.mu = permeability,
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
	local vars = table()

	for _,fluid in ipairs(fluids) do
	
		-- k is 0,1,2
		local function vorticity(k)
			local xs = {'x','y','z'}
			local i = (k+1)%3
			local j = (i+1)%3
			return {[fluid..' vorticity '..xs[k+1]] = template([[
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
]], 		{
				i = i,
				j = j,
				eqn = self,
				fluid = fluid,
			})}
		end
		
		vars:append{
			{[fluid..' rho'] = '*value = W.'..fluid..'_rho;'},
			{[fluid..' vx'] = '*value = W.'..fluid..'_v.x;'},
			{[fluid..' vy'] = '*value = W.'..fluid..'_v.y;'},
			{[fluid..' vz'] = '*value = W.'..fluid..'_v.z;'},
			{[fluid..' v'] = '*value = coordLen(W.'..fluid..'_v, x);'},
			{[fluid..' mx'] = '*value = U->'..fluid..'_m.x;'},
			{[fluid..' my'] = '*value = U->'..fluid..'_m.y;'},
			{[fluid..' mz'] = '*value = U->'..fluid..'_m.z;'},
			{[fluid..' m'] = '*value = coordLen(U->'..fluid..'_m, x);'},
			{[fluid..' P'] = '*value = W.'..fluid..'_P;'},
			{[fluid..' eInt'] = '*value = calc_'..fluid..'_eInt(W);'},
			{[fluid..' eKin'] = '*value = calc_'..fluid..'_eKin(W, x);'},
			{[fluid..' ePot'] = '*value = U->'..fluid..'_ePot;'},
			{[fluid..' eTotal'] = '*value = U->'..fluid..'_ETotal / W.'..fluid..'_rho;'},
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
		{Ex = '*value = U->epsE.x / U->eps;'},
		{Ey = '*value = U->epsE.y / U->eps;'},
		{Ez = '*value = U->epsE.z / U->eps;'},
		{E = '*value = sqrt(ESq(*U, x));'},
		{Bx = '*value = U->B.x;'},
		{By = '*value = U->B.y;'},
		{Bz = '*value = U->B.z;'},
		{B = '*value = sqrt(BSq(*U, x));'},
		{['EM energy'] = [[
	//*value = .5 * (coordLen(U->epsE) + coordLen(U->B) / U->mu);
	*value = .5 * (real3_len(U->epsE) + real3_len(U->B) / U->mu);
]]},
	}:append(table{'E','B'}:map(function(var,i)
		local field = assert( ({E='epsE', B='B'})[var] )
		return {['div '..var] = template([[
	*value = .5 * (0.
<?
for j=0,solver.dim-1 do
?>		+ (U[stepsize.s<?=j?>].<?=field?>.s<?=j?> 
			- U[-stepsize.s<?=j?>].<?=field?>.s<?=j?>
		) / grid_dx<?=j?>
<?
end 
?>	)<? 
if field == 'epsE' then 
?> / U->eps<?
end
?>;
]], {solver=self.solver, field=field})}
	end)):append{
		{BPot = '*value = U->BPot;'},
		{sigma = '*value = U->sigma;'},
		{eps = '*value = U->eps;'},
		{mu = '*value = U->mu;'},
	}

	return vars
end

-- can it be zero sized?
function TwoFluidEMHD:getEigenTypeCode()
	return template([[
typedef struct {
<? for _,fluid in ipairs(fluids) do ?>	
	// Roe-averaged vars
	real <?=fluid?>_rho;
	real3 <?=fluid?>_v;
	real <?=fluid?>_hTotal;

	// derived vars
	real <?=fluid?>_vSq;
	real <?=fluid?>_Cs;
<? end ?>
	real eps, mu; 
} <?=eqn.eigen_t?>;
]], {
	eqn = self,
	fluids = fluids,
})
end

function TwoFluidEMHD:getEigenDisplayVars()
	local vars = table()
	for _,fluid in ipairs(fluids) do
		vars:append{
			{[fluid..' rho'] = '*value = eigen->'..fluid..'_rho;'},
			{[fluid..' vx'] = '*value = eigen->'..fluid..'_v.x;'},
			{[fluid..' vy'] = '*value = eigen->'..fluid..'_v.y;'},
			{[fluid..' vz'] = '*value = eigen->'..fluid..'_v.z;'},
			{[fluid..' v'] = '*value = coordLen(eigen->'..fluid..'_v, xInt[0]);'},
			{[fluid..' hTotal'] = '*value = eigen->'..fluid..'_hTotal;'},
			{[fluid..' vSq'] = '*value = eigen->'..fluid..'_vSq;'},
			{[fluid..' Cs'] = '*value = eigen->'..fluid..'_Cs;'},
		}
	end
	vars:append{
		{eps = '*value = eigen->eps;'},
		{mu = '*value = eigen->mu;'},
	}
	return vars
end

return TwoFluidEMHD

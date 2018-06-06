--[[
maxwell but extended to include ADM metric influence
based on 2009 Alcubierre et al charged black holes
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local clnumber = require 'cl.obj.number'
local template = require 'template'

local GRMaxwell = class(Equation)
GRMaxwell.name = 'GRMaxwell'
GRMaxwell.numStates = 10
GRMaxwell.numWaves = 6
GRMaxwell.numIntStates = 6

GRMaxwell.consVars = {
	-- the vectors are contravariant with ^t component that are zero
	{epsE = 'real3'},
	{B = 'real3'},
	{BPot = 'real'},	-- used to calculate the B potential & remove div
	
	-- these aren't dynamic at all, but I don't want to allocate a separate buffer
	{sigma = 'real'},
	{eps = 'real'},
	{mu = 'real'},
}

GRMaxwell.mirrorVars = {{'epsE.x', 'B.x'}, {'epsE.y', 'B.y'}, {'epsE.z', 'B.z'}}

GRMaxwell.hasEigenCode = true
GRMaxwell.useSourceTerm = true
GRMaxwell.roeUseFluxFromCons = true

GRMaxwell.initStates = require 'init.euler'

function GRMaxwell:init(args)
	GRMaxwell.super.init(self, args)

	local NoDiv = require 'op.nodiv'
	self.solver.ops:insert(NoDiv{solver=self.solver})
end

function GRMaxwell:getCommonFuncCode()
	return template([[
real ESq(<?=eqn.cons_t?> U, sym3 gamma) { 
	return real3_weightedLenSq(U.epsE, gamma) / (U.eps * U.eps);
}

real BSq(<?=eqn.cons_t?> U, sym3 gamma) {
	return real3_weightedLenSq(U.B, gamma);
}

inline <?=eqn.prim_t?> primFromCons(<?=eqn.cons_t?> U, real3 x) { return U; }
inline <?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> W, real3 x) { return W; }
]], {
		eqn = self,
	})
end

GRMaxwell.initStateCode = [[
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
	global <?=eqn.cons_t?>* U = UBuf + index;

	//used
	real3 E = _real3(0,0,0);
	real3 B = _real3(0,0,0);
	real conductivity = 1.;
	
	//natural units say eps0 = 1/4pi, mu0 = 4pi
	//but waves don't make it to the opposite side...
	//mu0 eps0 = 1/c^2
	real permittivity = 1.; //1. / (4. * M_PI);
	real permeability = 1.; //4. * M_PI;
	
	//throw-away
	real rho = 0;
	real3 v = _real3(0,0,0);
	real P = 0;
	real ePot = 0;
	
	<?=code?>
	
	U->epsE = real3_scale(E, permittivity);
	U->B = B;
	U->BPot = 0;
	U->sigma = conductivity;
	U->eps = permittivity;
	U->mu = permeability;
}
]]

GRMaxwell.solverCodeFile = 'eqn/gr-maxwell.cl'

function GRMaxwell:getCalcEigenBasisCode() end

function GRMaxwell:getDisplayVars()
	local solver = self.solver
	return GRMaxwell.super.getDisplayVars(self):append{ 
		{E_u = '*valuevec = real3_scale(U->epsE, 1. / U->eps);', type='real3'},
	
		-- eps_ijk E^j B^k
		{S_l = '*valuevec = real3_scale(real3_cross(U->epsE, U->B), 1. / U->eps);', type='real3'},
		
		{energy = template([[
	<?=solver:getADMVarCode()?>
	*value = .5 * (real3_weightedLenSq(U->epsE, gamma) + real3_lenSq(U->B, gamma) / (U->mu * U->mu));
]], {solver=solver})},

	}
	--[=[ div E and div B ... TODO redo this with metric (gamma) influence 
	:append(table{'E','B'}:map(function(var,i)
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
	end))
	--]=]
end

GRMaxwell.eigenVars = table{
	{eps = 'real'},
	{mu = 'real'},
	{lambda = 'real'},
}

function GRMaxwell:eigenWaveCode(side, eig, x, waveIndex)
	if waveIndex == 0 or waveIndex == 1 then
		return '-'..eig..'->lambda'
	elseif waveIndex == 2 or waveIndex == 3 then
		return '0'
	elseif waveIndex == 4 or waveIndex == 5 then
		return eig..'->lambda'
	else
		error'got a bad waveIndex'
	end
end

function GRMaxwell:getFluxFromConsCode()
	return template([[
<? for side=0,solver.dim-1 do ?>
<?=eqn.cons_t?> fluxFromCons_<?=side?>(
	<?=eqn.cons_t?> U,
	real3 x<?=
	solver:getADMArgs()?>
) {
	<?=solver:getADMVarCode()?>
	
	real3 B_u = U.B;
	real3 epsE_u = U.epsE;

	//this only works for fluxFromCons
	real3 B_l = sym3_real3_mul(gamma, B_u);
	real3 epsE_l = sym3_real3_mul(gamma, epsE_u);
	
	real mu = U.mu;
	real eps = U.eps;
	
	return (<?=eqn.cons_t?>){
	<? if side == 0 then ?>
		.epsE = _real3(0., B_l.z / mu, -B_l.y / mu),
		.B = _real3(0., -epsE_l.z / eps, epsE_l.y / eps),
	<? elseif side == 1 then ?>
		.epsE = _real3(-B_l.z / mu, 0., B_l.x / mu),
		.B = _real3(epsE_l.z / eps, 0., -epsE_l.x / eps),
	<? elseif side == 2 then ?>
		.epsE = _real3(B_l.y / mu, -B_l.x / mu, 0.),
		.B = _real3(-epsE_l.y / eps, epsE_l.x / eps, 0.),
	<? end ?>
		.BPot = 0.,
		.sigma = 0.,
		.eps = 0.,
		.mu = 0.,
	};
}
<? end ?>
]], {
		eqn = self,
		solver = self.solver,
	})
end

return GRMaxwell

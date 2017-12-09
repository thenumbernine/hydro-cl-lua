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
	{epsE = 'real3'},
	{B = 'real3'},
	{BPot = 'real'},	-- used to calculate the B potential & remove div
	
-- there's really no reason to store these, they're not dynamic at all
	{sigma = 'real'},
	{eps = 'real'},
	{mu = 'real'},
}

GRMaxwell.mirrorVars = {{'epsE.x', 'B.x'}, {'epsE.y', 'B.y'}, {'epsE.z', 'B.z'}}

GRMaxwell.hasEigenCode = true
GRMaxwell.useSourceTerm = true
GRMaxwell.hasFluxFromCons = true

GRMaxwell.initStates = require 'init.euler'

function GRMaxwell:init(solver)
	GRMaxwell.super.init(self, solver)

	local NoDiv = require 'solver.nodiv'
	solver.ops:insert(NoDiv{solver=solver})
end

function GRMaxwell:getCodePrefix()
	return table{
		GRMaxwell.super.getCodePrefix(self),
		template([[
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
}),
	}:concat'\n'
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

function GRMaxwell:getSolverCode()
	return template(file['eqn/gr-maxwell.cl'], {
		eqn=self,
		solver=self.solver,
	})
end

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
	{alpha = 'real'},
	{det_gamma = 'real'},
	{detg_gUjj = 'real'},	-- g g^jj
}

return GRMaxwell

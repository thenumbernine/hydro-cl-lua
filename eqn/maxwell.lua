--[[
based on Trangenstein
but there's a more advanced version in the slides by Chatterjee "Finite Volume Time Domain (FVTD) Computations for Electromagnetic Scattering"
B,t + curl E = -rho H
D,t - curl H = -J_i - sigma E
D = eps' E
B = mu' H
sigma = omega eps"
rho = omega mu"
eps_r = eps' - j eps"
mu_r = mu' - j mu"
what is rho?
what is sigma? conductance
what is eps', eps", mu', mu", omega, j?
then using Shokin, Fedurok, Lebedev, Chubarov "Parallel FVTD for Solving Maxwell Equations in Dielectric-Metal Composite Media"
I get 
D,t - curl H = 0
B,t + curl E = 0
D = eps0 eps_r E
B = mu0 H
div D = 0
div B = 0
Drude model: eps(omega) = eps_inf - omega_p^2 / (omega (omega + i Gamma))
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local clnumber = require 'cl.obj.number'
local template = require 'template'

local Maxwell = class(Equation)
Maxwell.name = 'Maxwell'
Maxwell.numStates = 10
Maxwell.numWaves = 6
Maxwell.numIntStates = 6
Maxwell.mirrorVars = {{'epsE.x', 'B.x'}, {'epsE.y', 'B.y'}, {'epsE.z', 'B.z'}}

Maxwell.hasEigenCode = true
Maxwell.useSourceTerm = true
Maxwell.hasFluxFromCons = true

Maxwell.initStates = require 'init.euler'

function Maxwell:getTypeCode()
	return template([[
typedef union {
	real ptr[10];
	struct {
		real3 epsE;
		real3 B;
		real BPot;	//used to calculate the B potential & remove div
		real sigma;
		real eps;
		real mu;
	};
} <?=eqn.cons_t?>;

typedef <?=eqn.cons_t?> <?=eqn.prim_t?>;

]], {
	eqn = self,
})
end

function Maxwell:getCodePrefix()
	return table{
		Maxwell.super.getCodePrefix(self),
		template([[
//hmm, for E and B, even if the coord is 2D, we need all 3D components ...
//this means we need coordLen functions with guaranteed dimensions, including tangent spaces

real ESq(<?=eqn.cons_t?> U, real3 x) { 
	//return coordLenSq(U.epsE, x) / (U.eps * U.eps);
	return real3_lenSq(U.epsE) / (U.eps * U.eps);
}

real BSq(<?=eqn.cons_t?> U, real3 x) {
	//return coordLenSq(U.B, x);
	return real3_lenSq(U.B);
}

inline <?=eqn.prim_t?> primFromCons(<?=eqn.cons_t?> U, real3 x) { return U; }
inline <?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> W, real3 x) { return W; }
]], {
	eqn = self,
}),
	}:concat'\n'
end

function Maxwell:getInitStateCode()
	self.initState = self.initStates[self.solver.initStateIndex]
	assert(self.initState, "couldn't find initState "..self.solver.initStateIndex)	
	local code = self.initState.init(self.solver)	
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

]]..code..[[
	
	U->epsE = real3_scale(E, permittivity);
	U->B = B;
	U->BPot = 0;
	U->sigma = conductivity;
	U->eps = permittivity;
	U->mu = permeability;
}
]], {
	eqn = self,
})
end

function Maxwell:getSolverCode()
	return template(file['eqn/maxwell.cl'], {eqn=self, solver=self.solver})
end

function Maxwell:getDisplayVars()
	return table{
		{Ex = '*value = U->epsE.x / U->eps;'},
		{Ey = '*value = U->epsE.y / U->eps;'},
		{Ez = '*value = U->epsE.z / U->eps;'},
		{['|E|'] = '*value = sqrt(ESq(*U, x));'},
		{Bx = '*value = U->B.x;'},
		{By = '*value = U->B.y;'},
		{Bz = '*value = U->B.z;'},
		{['|B|'] = '*value = sqrt(BSq(*U, x));'},
		{energy = [[
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
end

function Maxwell:getVecDisplayVars()
	return table{
		{E = 'valuevec = real3_scale(U->epsE, 1. / U->eps);'},
		{B = 'valuevec = U->B;'},
		{S = 'valuevec = real3_scale(real3_cross(U->epsE, U->B), 1. / U->eps);'},
	}
end

-- can it be zero sized?
function Maxwell:getEigenTypeCode()
	return 'typedef struct { real eps, mu; } '..self.eigen_t..';'
end

function Maxwell:getEigenDisplayVars()
	return {
		{eps = '*value = eigen->eps;'},
		{mu = '*value = eigen->mu;'},
	}
end

return Maxwell

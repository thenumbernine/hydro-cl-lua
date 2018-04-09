--[[ 2000 Munz
TODO incorporate H and D fields
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'
local clnumber = require 'cl.obj.number'
local template = require 'template'

local Maxwell = class(Equation)
Maxwell.name = 'Maxwell'
Maxwell.numWaves = 8
Maxwell.numIntStates = 8

Maxwell.consVars = {
	{E = 'real3'},
	{B = 'real3'},
	{phi = 'real'},
	{psi = 'real'},
	{conductivity = 'real'},
	{charge = 'real'},
}

Maxwell.mirrorVars = {{'E.x', 'B.x'}, {'E.y', 'B.y'}, {'E.z', 'B.z'}}

Maxwell.hasEigenCode = true
Maxwell.useSourceTerm = true
Maxwell.hasFluxFromCons = true

Maxwell.initStates = require 'init.euler'

function Maxwell:init(solver)
	Maxwell.super.init(self, solver)
end

function Maxwell:getCodePrefix()
	return table{
		Maxwell.super.getCodePrefix(self),
		template([[
real ESq(<?=eqn.cons_t?> U, real3 x) { return real3_lenSq(U.E); }
real BSq(<?=eqn.cons_t?> U, real3 x) { return real3_lenSq(U.B); }
inline <?=eqn.prim_t?> primFromCons(<?=eqn.cons_t?> U, real3 x) { return U; }
inline <?=eqn.cons_t?> consFromPrim(<?=eqn.prim_t?> W, real3 x) { return W; }
]], {
	eqn = self,
}),
	}:concat'\n'
end

Maxwell.initStateCode = [[
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
	
	real permittivity = 1.;
	real permeability = 1.;
	
	//throw-away
	real rho = 0;
	real3 v = _real3(0,0,0);
	real P = 0;
	real ePot = 0;
	
	<?=code?>
	
	U->E = E;
	U->B = B;
	U->phi = 0;
	U->psi = 0;
	U->conductivity = conductivity;
	U->charge = 0;
}
]]

function Maxwell:getSolverCode()
	return template(file['eqn/glm-maxwell.cl'], {eqn=self, solver=self.solver})
end

function Maxwell:getDisplayVars()
	return Maxwell.super.getDisplayVars(self):append{ 
		{S = '*valuevec = real3_cross(U->E, U->B);', type='real3'},
		{energy = [[
	*value = .5 * (real3_lenSq(U->E) + real3_lenSq(U->B));
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
end

Maxwell.eigenVars = table{
	{nothing = 'real'},
}

function Maxwell:eigenWaveCodePrefix(side, eig, x, waveIndex)
	return ''
end

function Maxwell:eigenWaveCode(side, eig, x, waveIndex)
	return ({
		'-speedOfLight * chi',
		'-speedOfLight * gamma',
		'-speedOfLight',
		'-speedOfLight',
		'speedOfLight',
		'speedOfLight',
		'speedOfLight * gamma',
		'speedOfLight * chi',
	})[waveIndex+1] or error('got a bad waveIndex: '..waveIndex)
end

return Maxwell

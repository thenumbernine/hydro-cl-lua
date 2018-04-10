--[[ 2000 Munz
TODO incorporate H and D fields
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local range = require 'ext.range'
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
	{rhoCharge = 'real'},
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
	U->rhoCharge = 0;
}
]]

function Maxwell:getSolverCode()
	return template(file['eqn/glm-maxwell.cl'], {eqn=self, solver=self.solver})
end

-- k is 0,1,2
local function curl(eqn,k,result,field)
	local xs = {'x','y','z'}
	local i = (k+1)%3
	local j = (i+1)%3
	return {['curl '..field..' '..xs[k+1]] = template([[
	if (OOB(1,1)) {
		<?=result?> = 0.;
	} else {

<? if i+1 <= solver.dim then ?>
		global const <?=eqn.cons_t?>* Uim = U - stepsize.s<?=i?>;
		global const <?=eqn.cons_t?>* Uip = U + stepsize.s<?=i?>;
		real vim_j = Uim-><?=field?>.s<?=j?>;
		real vip_j = Uip-><?=field?>.s<?=j?>;
<? else ?>
		real vim_j = 0.;
		real vip_j = 0.;
<? end?>

<? if j+1 <= solver.dim then ?>
		global const <?=eqn.cons_t?>* Ujm = U - stepsize.s<?=j?>;
		global const <?=eqn.cons_t?>* Ujp = U + stepsize.s<?=j?>;
		real vjm_i = Ujm-><?=field?>.s<?=i?>;
		real vjp_i = Ujp-><?=field?>.s<?=i?>;
<? else ?>
		real vjm_i = 0.;
		real vjp_i = 0.;
<? end ?>

		<?=result?> = (vjp_i - vjm_i) / (2. * grid_dx<?=i?>)
				- (vip_j - vim_j) / (2. * grid_dx<?=j?>);
	}
]], {
		i = i,
		j = j,
		eqn = eqn,
		solver = eqn.solver,
		result = result,
		field = field,
	})}
end

function Maxwell:getDisplayVars()
	local vars = Maxwell.super.getDisplayVars(self):append{ 
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

	for _,field in ipairs{'E', 'B'} do
		local v = range(0,2):map(function(i) 
			return curl(self,i,'valuevec->s'..i,field) 
		end)
		vars:insert{['curl '..field]= template([[
	<? for i=0,2 do ?>{
		<?=select(2,next(v[i+1]))?>
	}<? end ?>
]], {v=v}), type='real3'}
	end

	return vars
end

Maxwell.eigenVars = table{
	{nothing = 'real'},
}

function Maxwell:eigenWaveCodePrefix(side, eig, x, waveIndex)
	return ''
end

function Maxwell:eigenWaveCode(side, eig, x, waveIndex)
	return ({
		'-speedOfLight * divPhiWavespeed',
		'-speedOfLight * divPsiWavespeed',
		'-speedOfLight',
		'-speedOfLight',
		'speedOfLight',
		'speedOfLight',
		'speedOfLight * divPsiWavespeed',
		'speedOfLight * divPhiWavespeed',
	})[waveIndex+1] or error('got a bad waveIndex: '..waveIndex)
end

return Maxwell

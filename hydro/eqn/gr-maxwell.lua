--[[
maxwell but extended to include ADM metric influence
based on 2009 Alcubierre et al charged black holes
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local Equation = require 'hydro.eqn.eqn'

local GRMaxwell = class(Equation)
GRMaxwell.name = 'GRMaxwell'
GRMaxwell.numStates = 10
GRMaxwell.numWaves = 6
GRMaxwell.numIntStates = 6

-- I'm working on making complex numbers exchangeable
GRMaxwell.scalar = 'real'
--GRMaxwell.scalar = 'cplx'	-- not supported in the least bit at the moment

GRMaxwell.vec3 = GRMaxwell.scalar..'3'
GRMaxwell.mat3x3 = GRMaxwell.scalar..'3x3'

GRMaxwell.susc_t = GRMaxwell.scalar

GRMaxwell.consVars = {
	-- the vectors are contravariant with ^t component that are zero
	{name='D', type=GRMaxwell.vec3},
	{name='B', type=GRMaxwell.vec3},
	{name='divBPot', type=GRMaxwell.scalar},	-- used to calculate the B potential & remove div
	
	-- these aren't dynamic at all, but I don't want to allocate a separate buffer
	{name='sigma', type=GRMaxwell.scalar},
	-- TODO: 1/eps, or 1/sqrt(eps) even better
	{name='eps', type=GRMaxwell.susc_t},
	{name='mu', type=GRMaxwell.susc_t},
}

GRMaxwell.roeUseFluxFromCons = true

GRMaxwell.initConds = require 'hydro.init.euler':getList()

function GRMaxwell:init(args)
	GRMaxwell.super.init(self, args)

	local NoDiv = require 'hydro.op.nodiv'()
	self.solver.ops:insert(NoDiv{solver=self.solver})
end

GRMaxwell.solverCodeFile = 'hydro/eqn/gr-maxwell.cl'

function GRMaxwell:getEnv()
	local env = GRMaxwell.super.getEnv(self)
	local scalar = self.scalar
	env.vec3 = self.vec3
	env.susc_t = self.susc_t
	env.scalar = scalar
	env.zero = scalar..'_zero'
	env.inv = scalar..'_inv'
	env.neg = scalar..'_neg'
	env.fromreal = scalar..'_from_real'
	env.add = scalar..'_add'
	env.sub = scalar..'_sub'
	env.mul = scalar..'_mul'
	env.real_mul = scalar..'_real_mul'
	return env
end

function GRMaxwell:getDisplayVars()
	local solver = self.solver
	return GRMaxwell.super.getDisplayVars(self):append{ 
		{name='E_u', code='value.vreal3 = real3_real_mul(U->D, 1. / U->eps);', type='real3'},
	
		-- eps_ijk E^j B^k
		{name='S_l', code='value.vreal3 = real3_real_mul(real3_cross(U->D, U->B), 1. / U->eps);', type='real3'},
		
		{name='energy', code=self:template[[
	<?=solver:getADMVarCode()?>
	value.vreal = .5 * (real3_weightedLenSq(U->D, gamma) + real3_lenSq(U->B, gamma) / (U->mu * U->mu));
]]},

	}
	--[=[ div E and div B ... TODO redo this with metric (gamma) influence 
	:append(table{'E','B'}:mapi(function(var)
		local field = assert( ({D='D', B='B'})[var] )
		return {name='div '..var, code=self:template([[
	value.vreal = .5 * (0.
<?
for j=0,solver.dim-1 do
?>		+ (U[solver->stepsize.s<?=j?>].<?=field?>.s<?=j?> 
			- U[-solver->stepsize.s<?=j?>].<?=field?>.s<?=j?>
		) / solver->grid_dx.s<?=j?>
<?
end 
?>	)<? 
if field == 'D' then 
?> / U->eps<?
end
?>;
]], {field=field})}
	end))
	--]=]
end

GRMaxwell.eigenVars = table{
	{name='eps', type='real'},
	{name='mu', type='real'},
	{name='lambda', type='real'},
}

function GRMaxwell:eigenWaveCode(args)
	local waveIndex = args.waveIndex
	if waveIndex == 0 or waveIndex == 1 then
		return '-('..args.eig..').lambda'
	elseif waveIndex == 2 or waveIndex == 3 then
		return '0'
	elseif waveIndex == 4 or waveIndex == 5 then
		return '('..args.eig..').lambda'
	else
		error'got a bad waveIndex'
	end
end

return GRMaxwell

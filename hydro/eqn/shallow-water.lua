local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'


local ShallowWater = class(Equation)
ShallowWater.name = 'ShallowWater'

ShallowWater.roeUseFluxFromCons = true

ShallowWater.initConds = require 'hydro.init.euler'


-- TODO primVars doesn't autogen displayVars, and therefore units doesn't matter
ShallowWater.primVars = {
	{name='h', type='real', units='m'},
	{name='v', type='real3', units='m/s', variance='u'},	-- contravariant
}

ShallowWater.consVars = {
	{name='h', type='real', units='m'},
	{name='m', type='real3', units='m^2/s', variance='u'},	-- contravariant
}

function ShallowWater:createInitState()
	ShallowWater.super.createInitState(self)
	self:addGuiVars{	
		{name='gravity', value=1, units='m/s^2'},
	}
end

function ShallowWater:initCodeModule_fluxFromCons()
	self.solver.modules:add{
		name = 'fluxFromCons',
		depends = {
			'solver.solver_t',
			'eqn.prim-cons',
			'coord.normal',
		},
		code = self:template[[
<?=eqn.cons_t?> fluxFromCons(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.cons_t?> U,
	real3 x,
	normalInfo_t n
) {
	<?=eqn.prim_t?> W = primFromCons(solver, U, x);
	real v_n = normalInfo_vecDotN1(n, W.v);
	real3 nU = normalInfo_u1(n);
	return (<?=eqn.cons_t?>){
		.h = U.h * v_n,
		.m = real3_add(
			real3_real_mul(U.m, v_n),	//h v^i v_n
			real3_real_mul(nU, .5 * solver->gravity * U.h * U.h)	//.5 g h^2 n^i
		),
	};
}
]],
	}
end

function ShallowWater:initCodeModuleCommon()
	self.solver.modules:add{
		name = 'eqn.common',
		depends = {
			'solver.solver_t',
			'eqn.cons_t',
		},
		code = self:template[[
real calc_C(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U) {
	return sqrt(solver->gravity * U.h);
}
]],
	}
end

function ShallowWater:initCodeModulePrimCons()
	self.solver.modules:add{
		name = 'eqn.prim-cons',
		depends = {
			'real3',
			'solver.solver_t',
			'eqn.prim_t',
			'eqn.cons_t',
		},
		code = self:template[[
<?=eqn.prim_t?> primFromCons(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U, real3 x) {
	return (<?=eqn.prim_t?>){
		.h = U.h,
		.v = real3_real_mul(U.m, 1. / U.h),
	};
}

<?=eqn.cons_t?> consFromPrim(constant <?=solver.solver_t?>* solver, <?=eqn.prim_t?> W, real3 x) {
	return (<?=eqn.cons_t?>){
		.h = W.h,
		.m = real3_real_mul(W.v, W.h),
	};
}
]],
	}

	-- only used by PLM
	self.solver.modules:add{
		name = 'eqn.dU-dW',
		depends = {
			'real3',
			'solver.solver_t',
			'eqn.prim_t',
			'eqn.cons_t',
		},
		code = self:template[[
<?=eqn.cons_t?> apply_dU_dW(
	constant <?=solver.solver_t?>* solver,
	<?=eqn.prim_t?> WA, 
	<?=eqn.prim_t?> W, 
	real3 x
) {
	return (<?=eqn.cons_t?>){
		.h = W.h,
		.m = real3_add(
			real3_real_mul(WA.v, W.h), 
			real3_real_mul(W.v, WA.h)),
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
		.h = U.h,
		.v = real3_sub(
			real3_real_mul(U.m, 1. / WA.h),
			real3_real_mul(WA.v, U.h / WA.h)),
	};
}
]],
	}
end

ShallowWater.initCondCode = [[
kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	
	// this is in all init/euler.lua
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
		.h = rho,
		.v = cartesianToCoord(v, x),
	};

	UBuf[index] = consFromPrim(solver, W, x);
}
]]

ShallowWater.solverCodeFile = 'hydro/eqn/shallow-water.cl'

function ShallowWater:getModuleDependsSolver() 
	return {
		'eqn.prim-cons',
		'eqn.common',
	}
end

ShallowWater.displayVarCodeUsesPrims = true

-- [=[
ShallowWater.predefinedDisplayVars = {
	'U h',
	'U v',
	'U v mag',
	'U v x',
	'U v y',
	'U v z',
}
--]=]

function ShallowWater:getDisplayVars()
	local vars = ShallowWater.super.getDisplayVars(self)
	
	vars:append{
		{name='v', code='value.vreal3 = W.v;', type='real3', units='m/s'},
		{name='wavespeed', code='value.vreal = calc_C(solver, *U);', units='m/s'},
	}

	vars:insert(self:createDivDisplayVar{
		field = 'v', 
		getField = function(U, j)
			return U..'->m.s'..j..' / '..U..'->h'
		end,
		units = '1/s',
	} or nil)

	vars:insert(self:createCurlDisplayVar{
		field = 'v',
		getField = function(U, j)
			return U..'->m.s'..j..' / '..U..'->h'
		end,
		units = '1/s',
	} or nil)

	return vars
end

ShallowWater.eigenVars = table{
	-- Roe-averaged vars
	{name='h', type='real', units='kg/m^3'},
	{name='v', type='real3', units='m/s'},
	-- derived vars
	{name='C', type='real', units='m/s'},
}

function ShallowWater:eigenWaveCodePrefix(n, eig, x)
	return self:template([[
	real C_nLen = <?=eig?>.C * normalInfo_len(<?=n?>);
	real v_n = normalInfo_vecDotN1(<?=n?>, <?=eig?>.v);
]], {
		eig = '('..eig..')',
		x = x,
		n = n,
	})
end

function ShallowWater:consWaveCodePrefix(n, U, x, W)
	return self:template([[
	real C_nLen = calc_C(solver, <?=U?>) * normalInfo_len(<?=n?>);
<? if not W then 
	W = 'W'
?>
	<?=eqn.prim_t?> W = primFromCons(solver, <?=U?>, <?=x?>);
<? end ?>
	real v_n = normalInfo_vecDotN1(n, <?=W?>.v);
]], {
		U = '('..U..')',
		W = W and '('..W..')' or nil,
		n = n,
		x = x,
	})
end

function ShallowWater:consWaveCode(n, U, x, waveIndex)
	if waveIndex == 0 then
		return '(v_n - C_nLen)'
	elseif waveIndex >= 1 and waveIndex <= 2 then
		return 'v_n'
	elseif waveIndex == 3 then
		return '(v_n + C_nLen)'
	end
	error'got a bad waveIndex'
end

ShallowWater.eigenWaveCode = ShallowWater.consWaveCode

return ShallowWater

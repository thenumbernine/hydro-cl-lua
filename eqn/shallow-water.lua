local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'
local xNames = require 'common'.xNames


local ShallowWater = class(Equation)
ShallowWater.name = 'ShallowWater'

ShallowWater.hasEigenCode = true
ShallowWater.hasFluxFromConsCode = true
ShallowWater.roeUseFluxFromCons = true

ShallowWater.initStates = require 'init.euler'


-- TODO primVars doesn't autogen displayVars, and therefore units doesn't matter
ShallowWater.primVars = {
	{name='h', type='real', units='m'},
	{name='v', type='real3', units='m/s', variance='u'},	-- contravariant
}

ShallowWater.consVars = {
	{name='h', type='real', units='m'},
	{name='m', type='real3', units='m/s', variance='u'},	-- contravariant
}

function ShallowWater:createInitState()
	ShallowWater.super.createInitState(self)
	self:addGuiVars{	
		{name='gravity', value=1, units='m/s^2'},
	}
end


function ShallowWater:getCommonFuncCode()
	return template([[
real calc_C(constant <?=solver.solver_t?>* solver, <?=eqn.cons_t?> U) {
	return sqrt(solver->gravity * U.h);
}
]], {
		solver = self.solver,
		eqn = self,
	})
end


function ShallowWater:getPrimConsCode()
	return template([[
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

]], {
		solver = self.solver,
		eqn = self,
	})
end

ShallowWater.initStateCode = [[
<? local xNames = require 'common'.xNames ?>
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	
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

ShallowWater.solverCodeFile = 'eqn/shallow-water.cl'

ShallowWater.displayVarCodeUsesPrims = true

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

		real vim_j = Uim->m.s<?=j?> / Uim->h;
		real vip_j = Uip->m.s<?=j?> / Uip->h;
		
		real vjm_i = Ujm->m.s<?=i?> / Ujm->h;
		real vjp_i = Ujp->m.s<?=i?> / Ujp->h;
		
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

	-- vorticity = [x ,y ,z] [v.x, v.y, v.z][
	-- = [v.z,y - v.y,z; v.x,z - v.z,x; v.y,x - v.x,y]
	if not require 'solver.meshsolver'.is(self.solver) then
		if self.solver.dim == 2 then
			vars:insert(vorticity(self,2,'value.vreal'))
		elseif self.solver.dim == 3 then
			local v = xNames:mapi(function(x,i)
				return vorticity(self,i-1,'value.vreal3.'..x) 
			end)
			vars:insert{name='vorticityVec', code=template([[
	<? for i=0,2 do ?>{
		<?=v[i+1].code?>
	}<? end ?>
]], {v=v}), type='real3', units='m/s^2'}
		end
	end

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
	return template([[
	real C_nLen = <?=eig?>.C * normalInfo_len(<?=n?>);
	real v_n = normalInfo_vecDotN1(<?=n?>, <?=eig?>.v);
]], {
		eig = '('..eig..')',
		x = x,
		n = n,
	})
end

function ShallowWater:consWaveCodePrefix(n, U, x, W)
	return template([[
	real C_nLen = calc_C(solver, <?=U?>) * normalInfo_len(<?=n?>);
<? if not W then 
	W = 'W'
?>
	<?=eqn.prim_t?> W = primFromCons(solver, <?=U?>, <?=x?>);
<? end ?>
	real v_n = normalInfo_vecDotN1(n, <?=W?>.v);
]], {
		eqn = self,
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

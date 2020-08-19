local class = require 'ext.class'
local table = require 'ext.table'
local EinsteinEqn = require 'hydro.eqn.einstein'
local makePartials = require 'hydro.eqn.makepartial'

local common = require 'hydro.common'
local xNames = common.xNames


local BSSNOKFiniteDifferenceEquationBase = class(EinsteinEqn)

-- seems all the hyperbolic formalisms listed in Alcubierre's book use alpha sqrt(gamma^ii) for the speed-of-light wavespeed
-- however the 2017 Ruchlin paper says to use gamma_ij
BSSNOKFiniteDifferenceEquationBase.cflMethod = '2008 Alcubierre'
--BSSNOKFiniteDifferenceEquationBase.cflMethod = '2013 Baumgarte et al, eqn 32'
--BSSNOKFiniteDifferenceEquationBase.cflMethod = '2017 Ruchlin et al, eqn 53'

function BSSNOKFiniteDifferenceEquationBase:init(args)
	self.cflMethod = args.cflMethod
	BSSNOKFiniteDifferenceEquationBase.super.init(self, args)
end

-- source: https://en.wikipedia.org/wiki/Finite_difference_coefficient 
-- derivCoeffs[derivative][order] = {coeffs...}
local derivCoeffs = {
	-- upwind 1st deriv coefficients
	{
		[2] = {[0] = -1, 1},
		[4] = {[0] = -3/2, 2, -1/2},
		[6] = {[0] = -11/6, 3, -3/2, 1/3},
		[8] = {[0] = -25/12, 4, -3, 4/3, -1/4},
	}
}

-- used by the BSSNOKFiniteDifference solvers
-- maybe a parent class for just them?
-- TODO higher orders plz
-- computes the partial between U[stepSize] and U[0]
function BSSNOKFiniteDifferenceEquationBase:makePartialUpwind(field, fieldType, nameOverride)
	local clnumber = require 'cl.obj.number'
	fieldType = fieldType or self:fieldTypeForVar(field)

	local suffix = 'l'
	if not field:find'_' then suffix = '_' .. suffix end

	local function add(x,y) return fieldType..'_add('..x..', '..y..')' end
	local function sub(x,y) return fieldType..'_sub('..x..', '..y..')' end
	local function real_mul(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
	local zero = fieldType..'_zero'
	local name = nameOverride or ('partial_'..field..suffix..'_upwind')
	local deriv = 1
	local order = 2
	local d1coeffs = assert(derivCoeffs[deriv][order], "couldn't find d/dx^"..(deriv/2).." coefficients of order "..order)
	
	local lines = table()
	if fieldType == 'real' then
		lines:insert('\treal3 '..name..';')
	elseif fieldType == 'cplx' then
		lines:insert('\tcplx3 '..name..';')
	elseif fieldType == 'real3' then
		lines:insert('\treal3x3 '..name..';')
	elseif fieldType == 'cplx3' then
		lines:insert('\tcplx3x3 '..name..';')
	else
		lines:insert('\t'..fieldType..' '..name..'[3];')
	end
	for i,xi in ipairs(xNames) do
		local namei
		if fieldType == 'real' 
		or fieldType == 'cplx'
		or fieldType == 'real3'
		or fieldType == 'cplx3'
		then
			namei = name..'.'..xi
		else
			namei = name..'['..(i-1)..']'
		end
		local expr = zero
		if i <= self.solver.dim then
			local U = 'U[0].'..field
			local expr = real_mul(U, clnumber(d1coeffs[0]))
			for j=1,#d1coeffs do
				local coeff = d1coeffs[j]
				local U = 'U['..j..' * updir.'..xi..' * solver->stepsize.'..xi..'].'..field
				expr = add(expr, real_mul(U, clnumber(coeff)))
			end
			expr = real_mul(expr, '1. / solver->grid_dx.'..xi)
		end
		expr = real_mul(expr, '(real)updir.'..xi)
		lines:insert('\t'..namei..' = '..expr..';')
	end
	return lines:concat'\n'
end

function BSSNOKFiniteDifferenceEquationBase:initCodeModule_setFlatSpace()
	local solver = self.solver
	solver.modules:add{
		name = 'mystery_C_U',
		code = [[
//TODO 2017 Ruchlin eqn. 8, what is C^i?
#define mystery_C_U	real3_zero
]],
	}

	solver.modules:add{
		name = 'setFlatSpace',
		depends = {'mystery_C_U'},
		code = self:template[[
void setFlatSpace(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* U,
	real3 x
) {
	U->alpha = 1.;
	U->beta_U = real3_zero;
	U->epsilon_LL = sym3_zero;
	U->W = 1;
	U->K = 0;
	U->ABar_LL = sym3_zero;

	//LambdaBar^i = Delta^i + C^i = Delta^i_jk gammaBar^jk = (connBar^i_jk - connHat^i_jk) gammaBar^jk + C^i
	//but when space is flat we have connBar^i_jk = connHat^i_jk and therefore Delta^i_jk = 0, Delta^i = 0, and LambdaBar^i = 0
	U->LambdaBar_U = mystery_C_U;

<? if eqn.useShift == 'HyperbolicGammaDriver' then
?>	U->B_U = real3_zero;
<? end
?>

	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0;
	U->M_U = real3_zero;
}
]],
	}
end

function BSSNOKFiniteDifferenceEquationBase:initCodeModuleCalcDT()
	self.solver.modules:add{
		name = 'eqn.calcDT',
		depends = {
			'eqn.common',
			'coord_sqrt_g_uu##',
		},
		code = self:template[[
kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 x = cellBuf[index].pos;
	const global <?=eqn.cons_t?>* U = UBuf + index;

<? if eqn.cflMethod == '2008 Alcubierre' then
?>	sym3 gamma_uu = calc_gamma_uu(U, x);
<? end 
?>
	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
<? 
if eqn.cflMethod == '2013 Baumgarte et al, eqn 32' then
	-- TODO if the grid is static then this only needs to be done once
	if side == 0 then 
?>		dt = (real)min(dt, solver->grid_dx.x);
<?	elseif side == 1 then 
?>		dt = (real)min(dt, .5 * solver->grid_dx.x * solver->grid_dx.y);
<? 	elseif side == 2 then 
?>		dt = (real)min(dt, .5 * solver->grid_dx.x * sin(.5 * solver->grid_dx.y) * solver->grid_dx.z);
<? 	end 
else
	if eqn.cflMethod == '2008 Alcubierre' then 
?>		//this is asserting alpha and W >0, which they should be
		real absLambdaMax = U->alpha * sqrt(gamma_uu.<?=sym(side+1,side+1)?>);
		real dx = solver->grid_dx.s<?=side?>;
		dt = (real)min(dt, dx / absLambdaMax);
<? 	elseif eqn.cflMethod == '2017 Ruchlin et al, eqn 53' then 
?>		// for wavespeeds alpha sqrt(gammaBar^ii)
		// and if we assume alpha > 1
		// and gamma_ii ~ 1/gamma^ii 
		// and gammaBar_ij ~ gammaHat_ij 
		// then the typical CFL equation: dt <= dx / lambda, lambda = alpha sqrt(gammaBar^ii)
		// turns into the SENR code: dt <= sqrt(gammaHat_ii) * dx
		real sqrt_gammaHat_ii = coord_sqrt_g_uu<?=side..side?>(x);
		real ds = sqrt_gammaHat_ii  * solver->grid_dx.s<?=side?>;
		dt = (real)min(dt, ds);
<? 	end
end
?>	}<? end ?>
	dtBuf[index] = dt;
}
]],
	}
end

function BSSNOKFiniteDifferenceEquationBase:getModuleDependsCommon()
	return {
		'setFlatSpace',
	}
end

function BSSNOKFiniteDifferenceEquationBase:getModuleDependsApplyInitCond()
	return {
		'setFlatSpace',
	}
end


return BSSNOKFiniteDifferenceEquationBase 

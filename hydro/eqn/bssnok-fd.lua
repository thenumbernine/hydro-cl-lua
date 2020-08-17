local class = require 'ext.class'
local table = require 'ext.table'
local EinsteinEqn = require 'hydro.eqn.einstein'
local makePartials = require 'hydro.eqn.makepartial'
local common = require 'hydro.common'
local template = require 'template'

local BSSNOKFiniteDifferenceEquationBase = class(EinsteinEqn)

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
	local xNames = common.xNames
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

function BSSNOKFiniteDifferenceEquationBase:initCodeModules()
	BSSNOKFiniteDifferenceEquationBase.super.initCodeModules(self)
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
		code = template([[
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
]], {eqn=self, solver=solver}),
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

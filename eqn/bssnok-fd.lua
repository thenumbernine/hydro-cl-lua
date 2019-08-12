local class = require 'ext.class'
local table = require 'ext.table'
local EinsteinEqn = require 'eqn.einstein'
local makePartials = require 'eqn.makepartial'
local common = require 'common'

local BSSNOKFiniteDifferenceEquationBase = class(EinsteinEqn)

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
	local d1coeffs = {1}	-- TODO more plz 
	
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
			for j,coeff in ipairs(d1coeffs) do
				local UR = 'Uup['..j..' * solver->stepsize.'..xi..'].'..field
				local UL = 'Uup['..(j-1)..' * solver->stepsize.'..xi..'].'..field
				expr = add(expr, real_mul(sub(UR, UL), clnumber(coeff)))
			end
			expr = real_mul(expr, '1. / solver->grid_dx.'..xi)
		end
		lines:insert('\t'..namei..' = '..expr..';')
	end
	return lines:concat'\n'
end

return BSSNOKFiniteDifferenceEquationBase 

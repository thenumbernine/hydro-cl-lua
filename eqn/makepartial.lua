--[[
TODO better think through the API
I do want 1st deriv functions that return vectors
I do want 2nd deriv functions that return degree-2's
but I would also like 4th deriv functions that return either vectors (one element per dim) 
--]]
local clnumber = require 'cl.obj.number'
local table = require 'ext.table'
local common = require 'hydro.common'

-- source: https://en.wikipedia.org/wiki/Finite_difference_coefficient 
-- derivCoeffs[derivative][order] = {coeffs...}
local derivCoeffs = {
	-- centered 1st deriv coefficients 
	{
		[2] = {.5},
		[4] = {2/3, -1/12},
		[6] = {3/4, -3/20, 1/60},
		[8] = {4/5, -1/5, 4/105, -1/280},
	},
	-- centered 2nd deriv coefficients
	{
		[2] = {[0] = -2, 1},
		[4] = {[0] = -5/2, 4/3, -1/12},
		[6] = {[0] = -49/18, 3/2, -3/20, 1/90},
		[8] = {[0] = -205/72, 8/5, -1/5, 8/315, -1/560},
	},
	-- centered 3rd deriv coefficients
	{
		[2] = {-1, 1/2},
		[4] = {-13/8, 1, -1/8},
		[6] = {-61/30, 169/120, -3/10, 7/240},
	},
	-- centered 4th deriv coefficients
	{
		[2] = {[0] = 6, -4, 1},
		[4] = {[0] = 28/3, -13/2, 2, -1/6},
		[6] = {[0] = 91/8, -122/15, 169/60, -2/5, 7/240},
	}
}

--[[
order = order of the finite difference
solver = solver
field = field of U buffer to use
fieldType = type of field
nameOverride = what name to use for the partial vars
	default: partial_field_l
--]]
local function makePartialRank1(deriv, order, solver, field, fieldType, nameOverride)
	local xNames = common.xNames
	local suffix = 'l'
	if not field:find'_' then suffix = '_' .. suffix end

	local function add(x,y) return fieldType..'_add('..x..', '..y..')' end
	local function sub(x,y) return fieldType..'_sub('..x..', '..y..')' end
	local function real_mul(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
	local zero = fieldType..'_zero'
	local name = nameOverride or ('partial_'..field..suffix)
	local d1coeffs = assert(derivCoeffs[deriv][order], "couldn't find d/dx^"..deriv.." coefficients of order "..order)
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
		if i <= solver.dim then
			for j,coeff in ipairs(d1coeffs) do
				local UR = 'U['..j..' * solver->stepsize.'..xi..'].'..field
				local UL = 'U[-'..j..' * solver->stepsize.'..xi..'].'..field
				expr = add(expr, real_mul(sub(UR, UL), clnumber(coeff)))
			end
			expr = real_mul(expr, '1. / solver->grid_dx.'..xi)
		end
		lines:insert('\t'..namei..' = '..expr..';')
	end
	return lines:concat'\n'
end

local function makePartialRank2(deriv, order, solver, field, fieldType, nameOverride)
	local xNames = common.xNames
	local symNames = common.symNames
	local from6to3x3 = common.from6to3x3
	local suffix = 'll'
	if not field:find'_' then suffix = '_' .. suffix end
	
	local function add(x,y) return fieldType..'_add('..x..', '..y..')' end
	local function sub(x,y) return fieldType..'_sub('..x..', '..y..')' end
	local function real_mul(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
	local zero = fieldType..'_zero'
	local name = nameOverride or ('partial2_'..field..suffix)
	local d1coeffs = assert(derivCoeffs[deriv/2][order], "couldn't find d/dx^"..(deriv/2).." coefficients of order "..order)
	local d2coeffs = assert(derivCoeffs[deriv][order], "couldn't find d/dx^"..deriv.." coefficients of order "..order)
	local lines = table()
	if fieldType == 'real' then
		lines:insert('\tsym3 '..name..';')
	else
		lines:insert('\t'..fieldType..' '..name..'[6];')
	end
	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi, xj = xNames[i], xNames[j]
		local nameij
		if fieldType == 'real' then
			nameij = name..'.'..xij
		else
			nameij = name..'['..(ij-1)..']'
		end
		if i > solver.dim or j > solver.dim then
			lines:insert('\t'..nameij..' = '..zero..';')
		elseif i == j then
			local expr = real_mul('U->'..field, d2coeffs[0])
			for k,coeff in ipairs(d2coeffs) do
				local UR = 'U['..k..' * solver->stepsize.s'..(i-1)..'].'..field
				local UL = 'U[-'..k..' * solver->stepsize.s'..(i-1)..'].'..field
				expr = add(expr, real_mul(add(UR, UL), clnumber(coeff)))
			end
			expr = real_mul(expr, '1. / (solver->grid_dx.'..xi..' * solver->grid_dx.'..xi..')')
			lines:insert('\t'..nameij..' = '..expr..';')
		else
			local expr = zero
			for k,coeff_k in ipairs(d1coeffs) do
				for l,coeff_l in ipairs(d1coeffs) do
					local URR = 'U['..k..' * solver->stepsize.'..xi..' + '..l..' * solver->stepsize.'..xj..'].'..field
					local ULL = 'U[-'..k..' * solver->stepsize.'..xi..' - '..l..' * solver->stepsize.'..xj..'].'..field
					local ULR = 'U[-'..k..' * solver->stepsize.'..xi..' + '..l..' * solver->stepsize.'..xi..'].'..field
					local URL = 'U['..k..' * solver->stepsize.'..xi..' - '..l..' * solver->stepsize.'..xi..'].'..field
					expr = add(expr, real_mul(sub(add(URR, ULL), add(ULR, URL)), clnumber(coeff_k * coeff_l)))
				end
			end
			expr = real_mul(expr, '1. / (solver->grid_dx.'..xi..' * solver->grid_dx.'..xj..')')
			lines:insert('\t'..nameij..' = '..expr..';')
		end
	end
	return lines:concat'\n'
end

local function makePartial1(...) return makePartialRank1(1, ...) end
local function makePartial2(...) return makePartialRank2(2, ...) end
local function makePartial3(...) return makePartialRank1(3, ...) end
local function makePartial4(...) return makePartialRank2(4, ...) end

return {
	derivCoeffs = derivCoeffs,
	makePartialRank1 = makePartialRank1,
	makePartialRank2 = makePartialRank2,
	makePartial1 = makePartial1,
	makePartial2 = makePartial2,
	makePartial3 = makePartial3,
	makePartial4 = makePartial4,
}

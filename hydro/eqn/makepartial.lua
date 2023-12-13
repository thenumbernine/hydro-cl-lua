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
-- separating out the denom is improving my numerical accuracy.  without doing so, bssnok-fd-num, cartesian, minkowski, RK4 diverges.
local derivCoeffs = {
	-- centered, antisymmetric 1st deriv coefficients 
	{
		[2] = {1, denom=2},
		[4] = {8, -1, denom=12},
		[6] = {45, -9, 1, denom=60},
		[8] = {672, -168, 32, -3, denom=840},
	},
	-- centered, symmetric 2nd deriv coefficients
	{
		[2] = {[0] = -2, 1},
		[4] = {[0] = -30, 16, -1, denom=12},
		[6] = {[0] = -490, 270, -27, 2, denom=180},
		[8] = {[0] = -205/72, 8/5, -1/5, 8/315, -1/560},
	},
	-- centered, antisymmetric 3rd deriv coefficients
	{
		[2] = {-2, 1, denom=2},
		[4] = {-13, 8, -1, denom=8},
		[6] = {-488, 338, -72, 7, denom=240},
	},
	-- centered, symmetric 4th deriv coefficients
	{
		[2] = {[0] = 6, -4, 1},
		[4] = {[0] = 56, -39, 12, -1, denom=6},
		[6] = {[0] = 2730, -1952, 676, -96, 7},
	}
}

--[[
order = order of the finite difference
solver = solver
field = field of U buffer to use
	optionally field can be function(offset) as a getter for the field at that offset of the U pointer
	if you use field as a function, be sure to provide nameOverride
fieldType = type of field
nameOverride = what name to use for the partial vars
	default: partial_field_l
srcName = which ptr to read offsets from, default is 'U'
--]]
local function makePartialRank1(deriv, order, solver, field, fieldType, nameOverride, srcName)
	local xNames = common.xNames
	local suffix = 'l'
	if type(field) == 'string' and not field:find'_' then
		suffix = '_' .. suffix
	end

	local add, sub, real_mul, zero
	if fieldType == 'real' then	--for readability, just inline the macros here
		add = function(x,y,xt,yt)
			return x..' + '..y, 'add'
		end
		sub = function(x,y,xt,yt)
			if yt == 'add' or yt == 'sub' then y = '('..y..')' end
			return x..' - '..y, 'sub'
		end
		real_mul = function(x,y,xt,yt)
			if xt == 'add' or xt == 'sub' then x = '('..x..')' end
			if yt == 'add' or yt == 'sub' then y = '('..y..')' end
			return x..' * '..y, 'mul'
		end
		zero = '0.'
	else
		add = function(x,y) return fieldType..'_add('..x..', '..y..')' end
		sub = function(x,y) return fieldType..'_sub('..x..', '..y..')' end
		real_mul = function(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
		zero = fieldType..'_zero'
	end

	srcName = srcName or 'U'
	local readField = 
		type(field) == 'string'
		and function(index) return srcName..'['..index..'].'..field end
		or field

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
		local expr, exprtype = zero, 'value'
		if i <= solver.dim then
			for j,coeff in ipairs(d1coeffs) do
				local UR = readField(j..' * solver->stepsize.'..xi)
				local UL = readField('-'..j..' * solver->stepsize.'..xi)
				
				local subexpr, subexprtype = sub(UR, UL, 'value', 'value')
				local mulexpr, mulexprtype = real_mul(subexpr, clnumber(coeff), subexprtype, 'value')
				expr, exprtype = add(expr, mulexpr, exprtype, mulexprtype)
			end
			local denom = 'solver->grid_dx.'..xi
			if d1coeffs.denom then
				denom = denom .. ' * ' .. clnumber(d1coeffs.denom)
			end
			expr, exprtype = real_mul(expr, '(1. / ('..denom..'))', exprtype, 'value')
		end
		lines:insert('\t'..namei..' = '..expr..';')
	end
	return lines:concat'\n'
end

local function makePartialRank2(deriv, order, solver, field, fieldType, nameOverride, srcName)
	local xNames = common.xNames
	local symNames = common.symNames
	local from6to3x3 = common.from6to3x3
	local suffix = 'll'
	if type(field) == 'string' and not field:find'_' then suffix = '_' .. suffix end
	
	local add, sub, real_mul, zero
	if fieldType == 'real' then	--for readability, just inline the macros here
		add = function(x,y,xt,yt)
			return x..' + '..y, 'add'
		end
		sub = function(x,y,xt,yt)
			if yt == 'add' or yt == 'sub' then y = '('..y..')' end
			return x..' - '..y, 'sub'
		end
		real_mul = function(x,y,xt,yt)
			if xt == 'add' or xt == 'sub' then x = '('..x..')' end
			if yt == 'add' or yt == 'sub' then y = '('..y..')' end
			return x..' * '..y, 'mul'
		end
		zero = '0.'
	else
		add = function(x,y) return fieldType..'_add('..x..', '..y..')' end
		sub = function(x,y) return fieldType..'_sub('..x..', '..y..')' end
		real_mul = function(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
		zero = fieldType..'_zero'
	end

	srcName = srcName or 'U'
	local readField = 
		type(field) == 'string'
		and function(index) return srcName..'['..index..'].'..field end
		or field

	local name = nameOverride or ('partial2_'..field..suffix)
	local d1coeffs = assert(derivCoeffs[deriv/2][order], "couldn't find d/dx^"..(deriv/2).." coefficients of order "..order)
	local d2coeffs = assert(derivCoeffs[deriv][order], "couldn't find d/dx^"..deriv.." coefficients of order "..order)
	local lines = table()
	if fieldType == 'real' then
		lines:insert('\treal3s3 '..name..';')
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
			local expr, exprtype = real_mul(readField'0', d2coeffs[0], 'value', 'value')
			for k,coeff in ipairs(d2coeffs) do
				local UR = readField(k..' * solver->stepsize.'..xi)
				local UL = readField('-'..k..' * solver->stepsize.'..xi)
				local addexpr, addexprtype = add(UR, UL, 'value', 'value')
				local mulexpr, mulexprtype = real_mul(addexpr, clnumber(coeff), addexprtype)
				expr, exprtype = add(expr, mulexpr, exprtype, mulexprtype)
			end
			local denom = 'solver->grid_dx.'..xi..' * solver->grid_dx.'..xi
			if d2coeffs.denom then
				denom = denom .. ' * ' .. clnumber(d2coeffs.denom)
			end
			expr, exprtype = real_mul(expr, '(1. / ('..denom..'))', exprtype, 'value')
			lines:insert('\t'..nameij..' = '..expr..';')
		else
			local expr, exprtype = zero, 'value'
			for k,coeff_k in ipairs(d1coeffs) do
				for l,coeff_l in ipairs(d1coeffs) do
					local URR = readField(k..' * solver->stepsize.'..xi..' + '..l..' * solver->stepsize.'..xj)
					local ULL = readField('-'..k..' * solver->stepsize.'..xi..' - '..l..' * solver->stepsize.'..xj)
					local ULR = readField('-'..k..' * solver->stepsize.'..xi..' + '..l..' * solver->stepsize.'..xj)
					local URL = readField(k..' * solver->stepsize.'..xi..' - '..l..' * solver->stepsize.'..xj)
					local addRRLLexpr, addRRLLtype = add(URR, ULL)
					local addLRRLexpr, addLRRLtype = add(ULR, URL)
					local subexpr, subexprtype = sub(addRRLLexpr, addLRRLexpr, addRRLLtype, addLRRLtype)
					local mulexpr, mulexprtype = real_mul(subexpr, clnumber(coeff_k * coeff_l), subexprtype, 'value')
					expr, exprtype = add(expr, mulexpr, exprtype, mulexprtype)
				end
			end
			local denom = 'solver->grid_dx.'..xi..' * solver->grid_dx.'..xj
			if d1coeffs.denom then
				denom = denom .. ' * ' .. clnumber(d1coeffs.denom)
			end
			expr, exprtype = real_mul(expr, '(1. / ('..denom..'))', exprtype, 'value')
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

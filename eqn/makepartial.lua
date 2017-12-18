local clnumber = require 'cl.obj.number'
local table = require 'ext.table'
require 'common'(_G)

local typeInfo = {
	real = {
		add = function(a,b) 
			if tonumber(a) == 0 then return b end
			if tonumber(b) == 0 then return a end
			return '('..a..') + ('..b..')' 
		end, 
		sub = function(a,b) return '('..a..') - ('..b..')' end, 
		scale = function(a,b) 
			if tonumber(a) == 1 then return b end
			if tonumber(b) == 1 then return a end
			return '('..a..') * ('..b..')' 
		end, 
		zero = '0',
	},
	real3 = {
		add = function(a,b) return 'real3_add('..a..', '..b..')' end,
		sub = function(a,b) return 'real3_sub('..a..', '..b..')' end,
		scale = function(a,b) return 'real3_scale('..a..', '..b..')' end,
		zero = '_real3(0,0,0)',
	},
	sym3 = {
		add = function(a,b) return 'sym3_add('..a..', '..b..')' end,
		sub = function(a,b) return 'sym3_sub('..a..', '..b..')' end,
		scale = function(a,b) return 'sym3_scale('..a..', '..b..')' end,
		zero = '_sym3(0,0,0,0,0,0)',
	},
	_3sym3 = {
		add = function(a,b) return '_3sym3_add('..a..', '..b..')' end,
		sub = function(a,b) return '_3sym3_sub('..a..', '..b..')' end,
		scale = function(a,b) return '_3sym3_scale('..a..', '..b..')' end,
		zero = '(_3sym3){.s={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}}',
	},
}

-- derivCoeffs[derivative][accuracy] = {coeffs...}
local derivCoeffs = {
	-- antisymmetric coefficients 
	{
		[2] = {.5},
		[4] = {2/3, -1/12},
		[6] = {3/4, -3/20, 1/60},
		[8] = {4/5, -1/5, 4/105, -1/280},
	},
	-- symmetric
	{
		[2] = {[0] = -2, 1},
		[4] = {[0] = -5/2, 4/3, -1/12},
		[6] = {[0] = -49/18, 3/2, -3/20, 1/90},
		[8] = {[0] = -205/72, 8/5, -1/5, 8/315, -1/560},
	},
}

local function makePartial(order, solver, field, fieldType)
	local suffix = 'l'
	if not field:find'_' then suffix = '_' .. suffix end
	local name = 'partial_'..field..suffix
	local fieldTypeInfo = assert(typeInfo[fieldType], "failed to find typeInfo for "..fieldType)
	local add, sub, scale, zero = fieldTypeInfo.add, fieldTypeInfo.sub, fieldTypeInfo.scale, fieldTypeInfo.zero
	local d1coeffs = assert(derivCoeffs[1][order], "couldn't find 1st derivative coefficients of order "..order)
	local lines = table{'\t'..fieldType..' '..name..'[3];\n'}
	for i,xi in ipairs(xNames) do
		local namei = name..'['..(i-1)..']'
		local expr = zero
		if i <= solver.dim then
			for j,coeff in ipairs(d1coeffs) do
				expr = add(expr, scale(sub(
						'U['..j..' * stepsize.'..xi..'].'..field,
						'U[-'..j..' * stepsize.'..xi..'].'..field
					), clnumber(coeff)))
			end
			expr = scale(expr, '1. / grid_dx'..(i-1))
		end
		lines:insert('\t'..namei..' = '..expr..';')
	end
	return lines:concat'\n'
end

local function makePartial2(order, solver, field, fieldType, nameOverride)
	local suffix = 'll'
	if not field:find'_' then suffix = '_' .. suffix end
	local name = nameOverride or ('partial2_'..field..suffix)
	local fieldTypeInfo = assert(typeInfo[fieldType], "failed to find typeInfo for "..fieldType)
	local add, sub, scale, zero = fieldTypeInfo.add, fieldTypeInfo.sub, fieldTypeInfo.scale, fieldTypeInfo.zero
	local d1coeffs = assert(derivCoeffs[1][order], "couldn't find 1st derivative coefficients of order "..order)
	local d2coeffs = assert(derivCoeffs[2][order], "couldn't find 2nd derivative coefficients of order "..order)
	local lines = table()
	lines:insert('\t'..fieldType..' '..name..'[6];')
	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi, xj = xNames[i], xNames[j]
		local nameij = name..'['..(ij-1)..']'
		if i > solver.dim or j > solver.dim then
			lines:insert('\t'..nameij..' = '..zero..';')
		elseif i == j then
			local expr = scale('U->'..field, d2coeffs[0])
			for k,coeff in ipairs(d2coeffs) do
				expr = add(
					expr, 
					scale(
						add(
							'U['..k..' * stepsize.s'..(i-1)..'].'..field,
							'U[-'..k..' * stepsize.s'..(i-1)..'].'..field),
						clnumber(coeff)))
			end
			expr = scale(expr, '1. / (grid_dx'..(i-1)..' * grid_dx'..(i-1)..')')
			lines:insert('\t'..nameij..' = '..expr..';')
		else
			local expr = zero
			for k,coeff_k in ipairs(d1coeffs) do
				for l,coeff_l in ipairs(d1coeffs) do
					expr = add(expr, scale(
						sub(
							add(
								'U['..k..' * stepsize.'..xi..' + '..l..' * stepsize.'..xj..'].'..field,
								'U[-'..k..' * stepsize.'..xi..' - '..l..' * stepsize.'..xj..'].'..field),
							add(
								'U[-'..k..' * stepsize.'..xi..' + '..l..' * stepsize.'..xi..'].'..field,
								'U['..k..' * stepsize.'..xi..' - '..l..' * stepsize.'..xi..'].'..field)), 
						clnumber(coeff_k * coeff_l)))
				end
			end
			expr = scale(expr, '1. / (grid_dx'..(i-1)..' * grid_dx'..(i-1)..')')
			lines:insert('\t'..nameij..' = '..expr..';')
		end
	end
	return lines:concat'\n'
end

return {
	makePartial = makePartial,
	makePartial2 = makePartial2,
}

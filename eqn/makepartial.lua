local clnumber = require 'cl.obj.number'
local table = require 'ext.table'
require 'common'(_G)

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
	
	local function add(x,y) return fieldType..'_add('..x..', '..y..')' end
	local function sub(x,y) return fieldType..'_sub('..x..', '..y..')' end
	local function real_mul(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
	local zero = fieldType..'_zero'

	local d1coeffs = assert(derivCoeffs[1][order], "couldn't find 1st derivative coefficients of order "..order)
	local lines = table{'\t'..fieldType..' '..name..'[3];\n'}
	for i,xi in ipairs(xNames) do
		local namei = name..'['..(i-1)..']'
		local expr = zero
		if i <= solver.dim then
			for j,coeff in ipairs(d1coeffs) do
				expr = add(expr, real_mul(sub(
						'U['..j..' * solver->stepsize.'..xi..'].'..field,
						'U[-'..j..' * solver->stepsize.'..xi..'].'..field
					), clnumber(coeff)))
			end
			expr = real_mul(expr, '1. / solver->grid_dx.'..xi)
		end
		lines:insert('\t'..namei..' = '..expr..';')
	end
	return lines:concat'\n'
end

local function makePartial2(order, solver, field, fieldType, nameOverride)
	local suffix = 'll'
	if not field:find'_' then suffix = '_' .. suffix end
	local name = nameOverride or ('partial2_'..field..suffix)
	
	local function add(x,y) return fieldType..'_add('..x..', '..y..')' end
	local function sub(x,y) return fieldType..'_sub('..x..', '..y..')' end
	local function real_mul(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
	local zero = fieldType..'_zero'

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
			local expr = real_mul('U->'..field, d2coeffs[0])
			for k,coeff in ipairs(d2coeffs) do
				expr = add(
					expr, 
					real_mul(
						add(
							'U['..k..' * solver->stepsize.s'..(i-1)..'].'..field,
							'U[-'..k..' * solver->stepsize.s'..(i-1)..'].'..field),
						clnumber(coeff)))
			end
			expr = real_mul(expr, '1. / (solver->grid_dx.'..xi..' * solver->grid_dx.'..xi..')')
			lines:insert('\t'..nameij..' = '..expr..';')
		else
			local expr = zero
			for k,coeff_k in ipairs(d1coeffs) do
				for l,coeff_l in ipairs(d1coeffs) do
					expr = add(expr, real_mul(
						sub(
							add(
								'U['..k..' * solver->stepsize.'..xi..' + '..l..' * solver->stepsize.'..xj..'].'..field,
								'U[-'..k..' * solver->stepsize.'..xi..' - '..l..' * solver->stepsize.'..xj..'].'..field),
							add(
								'U[-'..k..' * solver->stepsize.'..xi..' + '..l..' * solver->stepsize.'..xi..'].'..field,
								'U['..k..' * solver->stepsize.'..xi..' - '..l..' * solver->stepsize.'..xi..'].'..field)), 
						clnumber(coeff_k * coeff_l)))
				end
			end
			expr = real_mul(expr, '1. / (solver->grid_dx.'..xi..' * solver->grid_dx.'..xj..')')
			lines:insert('\t'..nameij..' = '..expr..';')
		end
	end
	return lines:concat'\n'
end

return {
	makePartial = makePartial,
	makePartial2 = makePartial2,
}

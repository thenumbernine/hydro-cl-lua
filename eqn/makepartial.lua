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

-- all this rescale stuff is specific to the bssn solver
local string = require 'ext.string'
local function deduceRescale(name, fieldType)
	local parts = string.split(name, '_')
	-- if there's no _ in the variable name then it's probably a scalar
	assert((#parts == 1) == (fieldType == 'real'), "variable "..name.." is of type "..fieldType..", but I think it should be a real")
	if #parts == 1 then return false, parts end
	return parts[#parts], parts
end

--[[
order = order of the finite difference
solver = solver
field = field of U buffer to use
fieldType = type of field
nameOverride = what name to use for the partial vars
	default: partial_field_l
rescale = nil to not rescale, or the tensor suffix to use rescaling
	real3: suffix can be u or l
	sym3: suffix can be uu or ll
	_3sym3: suffix can be uuu or lll
	(see coord/coord.lua for more on rescaling)
	Notice, rescale functions need the coordinate.  I'm assuming 'x' for now.
	Right now I'm just rescaling before calculating the partial, but not after, so the partial index is in coordinate form.
--]]
local function makePartial(order, solver, field, fieldType, nameOverride, rescale)
	local suffix = 'l'
	if not field:find'_' then suffix = '_' .. suffix end

	local function add(x,y) return fieldType..'_add('..x..', '..y..')' end
	local function sub(x,y) return fieldType..'_sub('..x..', '..y..')' end
	local function real_mul(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
	local zero = fieldType..'_zero'
	
	local name = nameOverride or ('partial_'..field..suffix)

	local parts
	if rescale == true then 
		rescale, parts = deduceRescale(field, fieldType) 
		-- lower the last _... of the name
		if rescale and not nameOverride then
			local a, b = name:match'^(.*)_(.-)$'
			if a and b then
				b = b:lower()
				name = a..'_'..b
			end
		end
	end

	local d1coeffs = assert(derivCoeffs[1][order], "couldn't find 1st derivative coefficients of order "..order)
	local lines = table{'\t'..fieldType..' '..name..'[3];\n'}
	for i,xi in ipairs(xNames) do
		local namei = name..'['..(i-1)..']'
		local expr = zero
		if i <= solver.dim then
			for j,coeff in ipairs(d1coeffs) do
				local UR = 'U['..j..' * solver->stepsize.'..xi..'].'..field
				local UL = 'U[-'..j..' * solver->stepsize.'..xi..'].'..field
				if rescale then
					UL = fieldType..'_rescaleToCoord_'..rescale..'('..UL..', x)'
					UR = fieldType..'_rescaleToCoord_'..rescale..'('..UR..', x)'
				end
				expr = add(expr, real_mul(sub(UR, UL), clnumber(coeff)))
			end
			expr = real_mul(expr, '1. / solver->grid_dx.'..xi)
			--[[
			if rescale then
				expr = fieldType..'_rescaleFromCoord_'..rescale..'('..expr..', x)'
			end
			--]]
		end
		lines:insert('\t'..namei..' = '..expr..';')
	end
	return lines:concat'\n'
end

local function makePartial2(order, solver, field, fieldType, nameOverride, rescale)
	local suffix = 'll'
	if not field:find'_' then suffix = '_' .. suffix end
	
	local function add(x,y) return fieldType..'_add('..x..', '..y..')' end
	local function sub(x,y) return fieldType..'_sub('..x..', '..y..')' end
	local function real_mul(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
	local zero = fieldType..'_zero'
	
	local name = nameOverride or ('partial2_'..field..suffix)

	local parts
	if rescale == true then 
		rescale, parts = deduceRescale(field, fieldType) 
		-- lower the last _... of the name
		if rescale and not nameOverride then
			local a, b = name:match'^(.*)_(.-)$'
			if a and b then
				b = b:lower()
				name = a..'_'..b
			end	
		end
	end

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
				local UR = 'U['..k..' * solver->stepsize.s'..(i-1)..'].'..field
				local UL = 'U[-'..k..' * solver->stepsize.s'..(i-1)..'].'..field
				if rescale then
					UL = fieldType..'_rescaleToCoord_'..rescale..'('..UL..', x)'
					UR = fieldType..'_rescaleToCoord_'..rescale..'('..UR..', x)'
				end
				expr = add(expr, real_mul(add(UR, UL), clnumber(coeff)))
			end
			expr = real_mul(expr, '1. / (solver->grid_dx.'..xi..' * solver->grid_dx.'..xi..')')
			--[[
			if rescale then
				expr = fieldType..'_rescaleFromCoord_'..rescale..'('..expr..', x)'
			end
			--]]
			lines:insert('\t'..nameij..' = '..expr..';')
		else
			local expr = zero
			for k,coeff_k in ipairs(d1coeffs) do
				for l,coeff_l in ipairs(d1coeffs) do
					local URR = 'U['..k..' * solver->stepsize.'..xi..' + '..l..' * solver->stepsize.'..xj..'].'..field
					local ULL = 'U[-'..k..' * solver->stepsize.'..xi..' - '..l..' * solver->stepsize.'..xj..'].'..field
					local ULR = 'U[-'..k..' * solver->stepsize.'..xi..' + '..l..' * solver->stepsize.'..xi..'].'..field
					local URL = 'U['..k..' * solver->stepsize.'..xi..' - '..l..' * solver->stepsize.'..xi..'].'..field
					if rescale then
						ULL = fieldType..'_rescaleToCoord_'..rescale..'('..ULL..', x)'
						ULR = fieldType..'_rescaleToCoord_'..rescale..'('..ULR..', x)'
						URL = fieldType..'_rescaleToCoord_'..rescale..'('..URL..', x)'
						URR = fieldType..'_rescaleToCoord_'..rescale..'('..URR..', x)'
					end
					expr = add(expr, real_mul(sub(add(URR, ULL), add(ULR, URL)), clnumber(coeff_k * coeff_l)))
				end
			end
			expr = real_mul(expr, '1. / (solver->grid_dx.'..xi..' * solver->grid_dx.'..xj..')')
			--[[
			if rescale then
				expr = fieldType..'_rescaleFromCoord_'..rescale..'('..expr..', x)'
			end
			--]]
			lines:insert('\t'..nameij..' = '..expr..';')
		end
	end
	return lines:concat'\n'
end

return {
	makePartial = makePartial,
	makePartial2 = makePartial2,
}

local ffi = require 'ffi'
local table = require 'ext.table'

local function countReals(vars)
	local structSize = 0
	for _,var in ipairs(vars) do
		local vartype
		if type(var) == 'string' then
			vartype = 'real'
		elseif type(var) == 'table' then
			vartype = select(2, next(var))
		end
		structSize = structSize + ffi.sizeof(vartype)
	end
	local numReals = structSize / ffi.sizeof'real'
	return numReals
end

local function makeStruct(name, vars)
	local numReals = countReals(vars)

	local lines = table()
	lines:insert'typedef union {'
	lines:insert('	real ptr['..numReals..'];')
	lines:insert('	struct {')
	for _,var in ipairs(vars) do
		if type(var) == 'string' then
			lines:insert('		real '..var..';')
			vartype = 'real'
		elseif type(var) == 'table' then
			local vn, vt = next(var)
			lines:insert('		'..vt..' '..vn..';')
		end
	end
	lines:insert('	};')
	lines:insert('} '..name..';')
	return lines:concat'\n'
end

return {
	makeStruct = makeStruct,
	countReals = countReals,
}

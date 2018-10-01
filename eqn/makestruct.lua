-- TODO move this out of eqn/
local ffi = require 'ffi'
local table = require 'ext.table'

local function countScalars(vars, scalar)
	scalar = scalar or 'real'
	local structSize = 0
	for _,var in ipairs(vars) do
		local vartype
		if type(var) == 'string' then
			vartype = scalar
		elseif type(var) == 'table' then
			vartype = select(2, next(var))
			assert(vartype, "expected vartype for var "..require 'ext.tolua'(var))
		end
		structSize = structSize + ffi.sizeof(vartype)
	end
	local numScalars = structSize / ffi.sizeof(scalar)
	return numScalars
end

local function makeStruct(name, vars, scalar)
	scalar = scalar or 'real'
	local numScalars = countScalars(vars, scalar)

	local lines = table()
	lines:insert'typedef union {'
	lines:insert('	'..scalar..' ptr['..numScalars..'];')
	lines:insert('	struct {')
	for _,var in ipairs(vars) do
		if type(var) == 'string' then
			lines:insert('		'..scalar..' '..var..';')
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

-- static so that no names overlap across all equations
local allnames = {}
function uniqueName(name)
	--[[ I'm using the base name in my typedefs per-source file
	if not allnames[name] then
		allnames[name] = true
		return name
	end
	--]]
	for i=2,math.huge do
		local try = name..'_'..i
		if not allnames[try] then
			allnames[try] = true
			return try
		end
	end
end

local function safeFFICDef(code)
	xpcall(function()
		ffi.cdef(code)
	end, function(msg)
		print(require 'template.showcode'(code))
		io.stderr:write(msg..'\n'..debug.traceback())
		os.exit(1)
	end)
end

return {
	makeStruct = makeStruct,
	safeFFICDef = safeFFICDef,
	countScalars = countScalars,
	uniqueName = uniqueName,
}

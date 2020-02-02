-- TODO move this out of eqn/
local ffi = require 'ffi'
local table = require 'ext.table'

--[[
eqn consVars/primVars structure:
vars = {
	{name=name, type=type, units=units}, ...
}
--]]

local function countScalars(vars, scalar)
	scalar = scalar or 'real'
	local structSize = 0
	for _,var in ipairs(vars) do
		structSize = structSize + ffi.sizeof(var.type)
	end
	local numScalars = structSize / ffi.sizeof(scalar)
	return numScalars
end

local function makeStruct(name, vars, scalar, dontUnion)
	scalar = scalar or 'real'
	local numScalars = countScalars(vars, scalar)

	local lines = table()
	local tab
	if dontUnion then
		lines:insert'typedef struct {'
		tab = '\t'
	else
		lines:insert'typedef union {'
		lines:insert('	'..scalar..' ptr['..numScalars..'];')
		lines:insert('	struct {')
		tab = '\t\t'
	end	
	for _,var in ipairs(vars) do
		
		lines:insert(
			tab
			..var.type
			
			-- fixing 'half' and 'double' alignment in solver_t
			-- dontUnion is only used by solver_t
			-- and solver_t is the only one with this C/CL alignment problem
			..(dontUnion and ' __attribute__ ((packed))' or '')
			
			..' '..var.name..';')
	end
	if not dontUnion then
		lines:insert('	};')
	end
	lines:insert('} '..name..';')
	return lines:concat'\n'
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
}

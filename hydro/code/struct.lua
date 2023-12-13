local ffi = require 'ffi'

local Struct = {}

function Struct:countScalars(scalar)
	scalar = scalar or 'real'
	local structSize = 0
	for _,var in ipairs(self.vars) do
		local res, err = xpcall(function()
			structSize = structSize + ffi.sizeof(var.type)
		end, function(err)
			return 'ffi.sizeof('..var.type..') failed for var '..require 'ext.tolua'(var)..'\n'
				..tostring(err)..'\n'
				..debug.traceback()
		end)
		if not res then error(err) end
	end
	local numScalars = structSize / ffi.sizeof(scalar)
	return numScalars
end

return Struct

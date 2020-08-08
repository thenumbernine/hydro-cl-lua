local ffi = require 'ffi'
local showcode = require 'template.showcode'
return function(code)
	xpcall(function()
		ffi.cdef(code)
	end, function(msg)
		print(showcode(code))
		io.stderr:write(msg..'\n'..debug.traceback())
		os.exit(1)
	end)
end

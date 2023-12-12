local ffi = require 'ffi'
local showcode = require 'template.showcode'
return function(code)
	if cmdline.debugcdefs then
		print('***** BEGIN CDEF *****')
		print(debug.traceback())
		print(showcode(code))
		print('***** END CDEF *****')
	end
	assert(xpcall(function()
		ffi.cdef(code)
	end, function(msg)
		return showcode(code)..'\n'
			..msg..'\n'..debug.traceback()
	end))
end

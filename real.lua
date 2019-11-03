local ffi = require 'ffi'

local realptr = ffi.new'realparam[1]'

local function real(x)
	realptr[0] = x
	return realptr
end

return real

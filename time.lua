require 'ffi.c.sys.time'
local ffi = require 'ffi'

local gettimeofday_tv = ffi.new'struct timeval[1]'
local function getTime()
	local results = ffi.C.gettimeofday(gettimeofday_tv, nil)
	return tonumber(gettimeofday_tv[0].tv_sec) + tonumber(gettimeofday_tv[0].tv_usec) / 1000000
end

local function getn(...)
	local t = {...}
	t.n = select('#', ...)
	return t
end
local function time(name, cb)
	print(name..'...')
	local startTime = getTime()
	local result = getn(cb())
	local endTime = getTime()
	print('...done '..name..' ('..(endTime - startTime)..'s)')
	return table.unpack(result, 1, result.n)
end

return {time, getTime}

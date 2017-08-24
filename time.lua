require 'ffi.c.sys.time'
local ffi = require 'ffi'

local getTime
if ffi.os == 'Windows' then
	-- in linux this is the live time, or something other than the actual time
	getTime = os.clock
else
	local gettimeofday_tv = ffi.new'struct timeval[1]'
	function getTime()
		local results = ffi.C.gettimeofday(gettimeofday_tv, nil)
		return tonumber(gettimeofday_tv[0].tv_sec) + tonumber(gettimeofday_tv[0].tv_usec) / 1000000
	end
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

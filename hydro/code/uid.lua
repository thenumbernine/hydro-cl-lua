--[[
returns the string of the unique id that Lua uses (the pointer of the object).
does so by clearing the mt and getting the tostring()
that means if the mt itself has a __metatable field then things could go wrong
--]]
--[[ method #1 ... use the memory address of the object
local function getuid(obj)
	local mt = getmetatable(obj)
	setmetatable(obj, nil)
	local uid = assert(tostring(obj):match'table: 0x(.*)')
	setmetatable(obj, mt)
	return uid
end
--]]
-- [[ method #2 ... use a cached pseudo-random # (reproducible) for each object
-- [=[ gcc rand from https://github.com/gcc-mirror/gcc/blob/master/libgfortran/intrinsics/rand.c
local ffi = require 'ffi'
local randseed = ffi.new('uint32_t', 123459876)
local function nextrand()
	randseed = ffi.cast('uint32_t', randseed * 16807) % 2147483647
	return randseed 
end
--]=]
local uidForObj = {}
local function getuid(obj)
	local mt = getmetatable(obj)
	setmetatable(obj, nil)
	local uidkey = assert(tostring(obj):match'table: 0x(.*)')
	setmetatable(obj, mt)
	
	local uid = uidForObj[uidkey]
	if not uid then
		uid = tostring(tonumber(nextrand()), 16)
		uidForObj[uidkey] = uid
	end
	return uid
end
--]]
return getuid

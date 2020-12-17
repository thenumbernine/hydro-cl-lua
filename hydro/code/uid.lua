--[[
returns the string of the unique id that Lua uses (the pointer of the object).
does so by clearing the mt and getting the tostring()
that means if the mt itself has a __metatable field then things could go wrong
--]]
local function getuid(obj)
	local mt = getmetatable(obj)
	setmetatable(obj, nil)
	local uid = assert(tostring(obj):match'table: 0x(.*)')
	setmetatable(obj, mt)
	return uid
end
return getuid

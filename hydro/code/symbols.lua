--[[
helper function for ...
assigns a 'symbolPrefix' to the obj based on the obj's Lua ptr (in code.uid)
then assigns a 'symbols' table with the prefix assigned to each
--]]
return function(obj, fields)
	local uid = require 'hydro.code.uid'(obj)
	obj.symbolPrefix = obj.name..'_'..uid..'_'

	obj.symbols = {}
	for _,field in ipairs(fields) do
		obj.symbols[field] = obj.symbolPrefix..field
	end
end

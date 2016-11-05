-- add line numbers to code output
local string = require 'ext.string'
local function showcode(code)
	return string.split(string.trim(code),'\n'):map(function(l,i) return i..':'..l end):concat'\n'
end
return showcode

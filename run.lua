#!/usr/bin/env luajit

for i=1,#arg do
	if arg[i] == 'clcpu' then
		-- use my own CPU driver
		package.loaded['ffi.OpenCL'] = require 'cl-cpu'
		break
	end
end

--[[ debugging
local tolua = require 'ext.tolua'
debug.sethook(function(hooktype)
	local info = debug.getinfo(2, 'nSl')
	if info.source and #info.source > 10 then info.source = info.source:sub(1, 10)..'...' end 
	print('hook', hooktype, tolua(info))
end, 'cr')
--]]

require 'app'():run()

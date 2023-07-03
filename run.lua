#!/usr/bin/env luajit

local useTrace
for i=1,#arg do
	if arg[i] == 'trace' then
		useTrace = true
	end
	if arg[i] == 'nojit' then
		jit.off(true, true)
	end
end

if useTrace then
	-- debugging
	local tolua = require 'ext.tolua'
	debug.sethook(function(hooktype)
		local info = debug.getinfo(2, 'nSl')
		if info.source and #info.source > 10 then info.source = info.source:sub(1, 10)..'...' end
		print('hook', hooktype, tolua(info))
	end, 'cr')
end

return require 'hydro.app'():run()

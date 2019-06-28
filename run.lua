#!/usr/bin/env luajit

for i=1,#arg do
	if arg[i] == 'clcpu' then
		-- use my own CPU driver
		package.loaded['ffi.OpenCL'] = require 'cl-cpu'
		break
	end
end

require 'app'():run()

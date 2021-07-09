#!/usr/bin/env luajit
require 'ext'
local exec = require 'exec'	-- TODO put this somewhere everyone can get it

local ffi = require 'ffi'
local unistd = require 'ffi.c.unistd'	-- getcwd, chdir
require 'ffi.c.stdlib'	-- free
local ccwd = unistd.getcwd(nil, 0)
local cwd = ffi.string(ccwd)
ffi.C.free(ccwd)
assert(unistd.chdir'../..' == 0)
print('cwd = '..cwd)

-- for how long should this run?
-- TODO pick coord ... and pick vectorComponent
local cmd = 'luajit run.lua sys=console "config='..cwd..'\\config.lua" coord=cartesian exitTime=.5 "trackvars=U E_g mag" > "'..cwd..'\\out-cartesian.txt"'

exec(cmd)

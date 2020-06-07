#!/usr/bin/env luajit
local ffi = require 'ffi'
local unistd = require 'ffi.c.unistd'
require 'ffi.c.stdlib'
local dirp = unistd.getcwd(nil, 0)
local dir = ffi.string(dirp)
ffi.C.free(dirp)

-- chdir to the base
-- (another alternative would be to execute this script from the base)
-- I think that's what the Relativity project did
unistd.chdir'../..'

local table = require 'ext.table'
local class = require 'ext.class'

local cols = table()
local rows = table()

local f = io.open(dir..'/var-ranges.txt', 'w')

cmdline = {sys='console'}

local App = class(require 'hydro.app')

function App:setup()
	local solver = require 'hydro.solver.roe'{
		app = self, 
		eqn = 'euler',
		dim = 2,
		fluxLimiter = 'superbee',
		coord = 'cartesian',
		mins = {-1, -1, -1},
		maxs = {1, 1, 1},
		gridSize = {128, 128, 128},
		boundary = {
			xmin='mirror',
			xmax='mirror',
			ymin='mirror',
			ymax='mirror',
			zmin='mirror',
			zmax='mirror',
		},
		initState = 'Sod',
		integrator = 'forward Euler',
	}
	self.solvers:insert(solver)

	-- start the solver off
	self.running = true
	self.exitTime = 3

	self.trackVars = table()

	for _,var in ipairs(solver.displayVars) do
		if var.name:sub(1,2) == 'U ' 
		and not var.vectorField
		then
			self.trackVars:insert(var)
		end
	end

	-- initialize what variables to track
	f:write'#t'
	for _,var in ipairs(self.trackVars) do
		f:write('\t',var.name..'_min')
		f:write('\t',var.name..'_max')
		f:write('\t',var.name..'_avg')
	end
	f:write'\n'
	f:flush()
end

function App:update()
	App.super.update(self)

	-- also TODO - make dumpFile more modular
	-- with headers and getters
	local solver = self.solvers[1]
	f:write(solver.t)
	for _,var in ipairs(self.trackVars) do
		local ymin, ymax, yavg = solver:calcDisplayVarRangeAndAvg(var)
		f:write(('\t%.16e'):format(ymin))
		f:write(('\t%.16e'):format(ymax))
		f:write(('\t%.16e'):format(yavg))
	end
	f:write'\n'
	f:flush()
end

App():run()

f:close()

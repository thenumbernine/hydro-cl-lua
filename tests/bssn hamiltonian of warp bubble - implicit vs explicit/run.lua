#!/usr/bin/env luajit

local ffi = require 'ffi'
require 'ffi.c.unistd'
require 'ffi.c.stdlib'
local dirp = ffi.C.getcwd(nil, 0)
local dir = ffi.string(dirp)
ffi.C.free(dirp)

-- chdir to the base
-- (another alternative would be to execute this script from the base)
-- I think that's what the Relativity project did
ffi.C.chdir'../..'

local table = require 'ext.table'
local class = require 'ext.class'

__disableGUI__ = true -- set before require 'app'

local cols = table()
local rows = table()

for _,info in ipairs{
	{'fe', 'forward Euler'},
	{'be', 'backward Euler'},
} do
	local suffix, integrator= table.unpack(info)
	local f = io.open(dir..'/var-ranges-'..suffix..'.txt', 'w')

	local App = class(require 'app')
	
	function App:setup()
		local solver = require 'solver.z4c-fd'{
			app = self, 
			dim = 2,
			fluxLimiter = 'superbee',
			coord = 'cartesian',
			mins = {-1, -1, -1},
			maxs = {1, 1, 1},
			gridSize = {16, 16, 16},
			initState = 'Alcubierre warp bubble',
			integrator = integrator,
		}
		self.solvers:insert(solver)

		-- start the solver off
		self.running = true
		self.exitTime = .5

		
		-- track var system
		-- is currently intertwined with the display system

		self.trackVars = table()

		-- so this is interesting
		-- calcDisplayVarRange uses a kernel that is only created upon necessity
		-- which seems slick, except that 
		-- now I want the kernel for output and not for displaying
		-- but lazy me will enable it for displaying for the time being
		for _,var in ipairs(solver.displayVars) do
			if var.name:sub(1,2) == 'U ' 
			-- seems to be crashing when reducing vector fields ...
			and not var.vectorField
			then
				self.trackVars:insert(var)
			end
			var.enabled = false
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
end

print'done!'

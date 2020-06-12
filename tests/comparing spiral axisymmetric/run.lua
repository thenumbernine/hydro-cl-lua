#!/usr/bin/env luajit

local ffi = require 'ffi'
require 'ffi.c.stdlib'		-- free
local unistd = require 'ffi.c.unistd'
local class = require 'ext.class'

-- save the cwd and chdir to ../..
local rundirp = unistd.getcwd(nil, 0)
local rundir = ffi.string(rundirp)
ffi.C.free(rundirp)
local resultsDir = 'results'
os.execute('mkdir "'..rundir..'/'..resultsDir..'" 2> '..(ffi.os == 'Windows' and 'NIL' or '/dev/null'))
unistd.chdir'../..'

-- set this global to have hydro run in console mode
-- working on doing this cleaner...
cmdline = {
	sys = 'console',
	--showfps = true,
}

local HydroApp = class(require 'hydro.app')

function HydroApp:setup(args)
	local MeshSolver = require 'hydro.solver.meshsolver'

	function MeshSolver:update(...)
		MeshSolver.super.update(self, ...)
		
		local var = self.displayVarForName['U v']
		local component = self.displayComponentFlatList[var.component]
		assert(self:isVarTypeAVectorField(component.type))
		local vMin, vMax, vAvg = self:calcDisplayVarRangeAndAvg(var, component.magn)

		print(self.t, vMin, vAvg, vMax)
	end

	self.solvers:insert(
		MeshSolver{
			app = self,
			integrator = 'forward Euler',
			cfl = .5,
			eqn = 'euler',
			initState = 'spiral',
			-- no / donor cell flux limiter
			-- cartesian / holonomic vector components
			-- legacy to gridsolvers.  mesh doesn't need this:
			dim = 2,
			-- mesh-specific params:
			mesh = {
				type = 'polar2d',
				size = {64, 64},
			}
		}
	)
end

HydroApp():run()

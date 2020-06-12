#!/usr/bin/env luajit

local ffi = require 'ffi'
require 'ffi.c.stdlib'		-- free
local unistd = require 'ffi.c.unistd'
local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local io = require 'ext.io'
local tolua = require 'ext.tolua'

-- save the cwd and chdir to ../..
local rundirp = unistd.getcwd(nil, 0)
local rundir = ffi.string(rundirp)
ffi.C.free(rundirp)

--local resultsDir = 'results'
--os.execute('mkdir "'..rundir..'/'..resultsDir..'" 2> '..(ffi.os == 'Windows' and 'NIL' or '/dev/null'))

unistd.chdir'../..'

-- set this global to have hydro run in console mode
-- working on doing this cleaner...
cmdline = {
	--showfps = true,
	sys = 'console',
	exitTime = 10,
}

-- TODO don't build all at once
-- in fact, put the config stuff in a separate launcher
-- so the app's cdef namespace doesn't get overly cluttered
--
-- another TODO for hydro.app ... reuse names for matching ctypes
--
local data = table()

--local f = assert(io.open('v.txt', 'w'))

local HydroApp = class(require 'hydro.app')

function HydroApp:setup(args)
	
	local function addSolver(identifier, solver)
		local oldUpdate = solver.update
		
		function solver:update(...)
			oldUpdate(self, ...)
			
			local var = self.displayVarForName['U v']
			local component = self.displayComponentFlatList[var.component]
			assert(self:isVarTypeAVectorField(component.type))
			local vMin, vMax, vAvg = self:calcDisplayVarRangeAndAvg(var, component.magn)

			--f:write(self.t,'\t', vMin,'\t', vAvg,'\t', vMax,'\n')
			--f:flush()
			local datai = data[self.identifier]
			if not datai then
				datai = table()
				data[self.identifier] = datai
			end
			datai:insert{self.t, vMin, vAvg, vMax}
			
			local floorT = math.floor(self.t)
			if not lastTime or floorT ~= lastTime then
				lastTime = floorT
				print(self.identifier, self.t, vMin, vAvg, vMax)
			end
		end

		self.solvers:insert(solver)
		solver.identifier = identifier
	end

--[[ mesh with cartesian components
-- I'm using this as my golden standard for curvilinear grids
	addSolver('mesh', require 'hydro.solver.meshsolver'{
		app = self,
		integrator = 'forward Euler',
		cfl = .25,
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
		},
	})
--]]
--[[ grid with cartesian components
-- this is nearly identical to the above, which is good.	
	addSolver('grid-cartesian', require 'hydro.solver.roe'{
		app = self,
		integrator = 'forward Euler',
		cfl = .25,
		eqn = 'euler',
		initState = 'spiral',
		coord = 'cylinder',
		coordArgs = {vectorComponent = 'cartesian'},
		dim = 2,
		mins = {.1, 0, -.5},
		maxs = {1, 2*math.pi, .5},
		gridSize = {64, 64, 1},
		boundary = {
			xmin='freeflow',
			xmax='freeflow',
			ymin='periodic',
			ymax='periodic',
			zmin='freeflow',
			zmax='freeflow',
		},
	})
--]]
--[[ grid with orthonormal grid-aligend components
-- this starts to deviate ...
	addSolver('grid-orthonormal', require 'hydro.solver.roe'{
		app = self,
		integrator = 'forward Euler',
		cfl = .25,
		eqn = 'euler',
		initState = 'spiral',
		coord = 'cylinder',
		coordArgs = {vectorComponent = 'anholonomic'},
		dim = 2,
		mins = {.1, 0, -.5},
		maxs = {1, 2*math.pi, .5},
		gridSize = {64, 64, 1},
		boundary = {
			xmin='freeflow',
			xmax='freeflow',
			ymin='periodic',
			ymax='periodic',
			zmin='freeflow',
			zmax='freeflow',
		},
	})
--]]
--[[ grid with grid coordinate components
-- this encounters numericals errors
	addSolver('grid-coordinate', require 'hydro.solver.roe'{
		app = self,
		integrator = 'forward Euler',
		cfl = .25,
		eqn = 'euler',
		initState = 'spiral',
		coord = 'cylinder',
		coordArgs = {vectorComponent = 'holonomic'},
		dim = 2,
		mins = {.1, 0, -.5},
		maxs = {1, 2*math.pi, .5},
		gridSize = {64, 64, 1},
		boundary = {
			xmin='freeflow',
			xmax='freeflow',
			ymin='periodic',
			ymax='periodic',
			zmin='freeflow',
			zmax='freeflow',
		},
	})
--]]
end

function HydroApp:requestExit()
	unistd.chdir(rundir)

	for identifier,datai in pairs(data) do
		file['results-'..identifier..'.txt'] = 
			'#t	mesh-v-min	mesh-v-avg	mesh-v-max\n'
			..datai:mapi(function(l) return table.concat(l, '\t') end):concat'\n'
			..'\n'
	end

	--os.execute('gnuplot plot.gnuplot')
	dofile'plot.lua'
	
	HydroApp.super.requestExit(self)
end

HydroApp():run()

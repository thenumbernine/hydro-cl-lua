#!/usr/bin/env luajit
require 'ext'
local ffi = require 'ffi'
local unistd = require 'ffi.c.unistd'
require 'ffi.c.stdlib'
local dirp = unistd.getcwd(nil, 0)
local dir = ffi.string(dirp)
ffi.C.free(dirp)
unistd.chdir'../..'

-- TODO make App accept this as a ctor argument
cmdline = {sys='console'}

local App = class(require 'hydro.app')

function App:setup()
	local solver = require 'solver.roe'{
		eqn = 'bssnok-fd',
		app = self, 
		dim = 3,
		integrator = 'backward Euler',
		fluxLimiter = 'superbee',
		coord = 'cartesian',
		mins = {-1, -1, -1},
		maxs = {1, 1, 1},
		gridSize = {32, 32, 32},
		boundary = {
			xmin='mirror',
			xmax='mirror',
			ymin='mirror',
			ymax='mirror',
			zmin='mirror',
			zmax='mirror',
		},
		initState = 'Alcubierre warp bubble',
	}

	self.solvers:insert(solver)

	self.thread = (function()
		self.running = true

		solver:save(dir..'/r1f1')
		solver:update()
		coroutine.yield()
		solver:save(dir..'/r1f2')

		-- now run until the dt is nan

		--while math.isfinite(solver.dt) do
		while solver.dt < 100 do
			solver:update()
			print(require 'ext.tolua'{t=solver.t, dt=solver.dt})
			coroutine.yield()
		end

		print'resetting...'
		solver:resetState()
		coroutine.yield()
		
		solver:save(dir..'/r2f1')
		solver:update()
		coroutine.yield()
		solver:save(dir..'/r2f2')

		self:requestExit()
	end):wrap()

-- [[ no gui
	while not self.done do self.thread() end
	os.exit()
--]]
end

--[[ with gui
function App:update()
	App.super.update(self)
	self.thread() 
end
--]]

App():run()

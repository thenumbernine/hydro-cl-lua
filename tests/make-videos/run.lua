#!/usr/bin/env luajit
--[[
this will make movies of predefined configurations
it expects ffmpeg to be installed
--]]

local ffi = require 'ffi'
package.cpath = package.cpath .. [[;C:\Users\moorece\lib\luajit-2.1.0-beta3\?.dll]]
local lfs = require 'lfs'
local rundir = lfs.currentdir()
lfs.chdir'../..'

local sdl = require 'ffi.sdl'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'

for _,setup in ipairs{
	function(self, args)
		local solverName = 'roe'
		local args = {
			app = self,
			eqn = 'euler',
			dim = 1,
			integrator = 'forward Euler',
			fluxLimiter = 'superbee',
			geometry = 'cartesian',
			mins = {-1,-1,-1},
			maxs = {1,1,1},
			gridSize = {256,1,1},
			boundary = {
				xmin = 'freeflow',
				xmax = 'freeflow',
				ymin = 'freeflow',
				ymax = 'freeflow',
				zmin = 'freeflow',
				zmax = 'freeflow',
			},
			initState = 'Sod',
		}
		self.solvers:insert(require('solver.'..solverName)(args))
		
		self.destMovieName = table{
			'eqn='..args.eqn,
			'solver='..solverName,
			'integrator='..args.integrator,
		
			-- TODO this only if usePLM is false 
			'fluxLimiter='..args.fluxLimiter,
			
			'init='..args.initState,
			'gridSize='..range(args.dim):map(function(i) return args.gridSize[i] end):concat'x',
		}:concat', '..'.mp4'
		print(self.destMovieName)
	end,
} do
	local App = class(require 'app')
	function App:setup(args)
		sdl.SDL_SetWindowSize(self.window, 1280, 720)
		setup(self, args)
		self.running = true
		self.exitTime = .1
		self.createAnimation = true
	end
	local app  = App()
	app:run()

	local function run(...)
		print(...)
		return os.execute(...)
	end

	-- once we're done ...
	local ext = app.screenshotExts[app.screenshotExtIndex]
	local dir = 'screenshots/'..app.screenshotDir
	assert(run('ffmpeg -y -i "'..dir..'/%05d.'..ext..'" "'..rundir..'/'..app.destMovieName..'"'))
	for i=0,app.screenshotIndex-1 do
		os.remove(dir..('/%05d.'):format(i)..ext)
	end
	local sep = ffi.os == 'Windows' and '\\' or '/'
	run('rmdir '..dir:gsub('/', sep))
end

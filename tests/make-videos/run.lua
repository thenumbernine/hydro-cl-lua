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
local io = require 'ext.io'

__disableGUI__ = true	-- set this before require 'app'

local function run(...)
	print(...)
	return os.execute(...)
end

local configurations = {
	{eqn='euler', gridSize={256}, initState='Sod', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
	{eqn='euler', gridSize={256}, initState='Sod', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
	{eqn='euler', gridSize={256}, initState='Sedov', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
	{eqn='euler', gridSize={256}, initState='Sedov', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
	{eqn='euler', gridSize={256}, initState='Sedov', solver='roe', integrator='backward Euler', usePLM='plm-eig-prim-ref'},		-- b.e. gets nans
	{eqn='euler', gridSize={256}, initState='Sedov', solver='hll', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
	{eqn='euler', gridSize={256}, initState='self-gravitation test 1', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
	{eqn='euler', gridSize={256}, initState='self-gravitation test 1', solver='hll', integrator='forward Euler', fluxLimiter='superbee'},

	{eqn='euler', gridSize={256,256}, initState='Sod', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
	{eqn='euler', gridSize={256,256}, initState='Sod', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
	{eqn='euler', gridSize={256,256}, initState='Sod', solver='hll', integrator='forward Euler', fluxLimiter='superbee'},
	{eqn='euler', gridSize={256,256}, initState='Sod', solver='hll', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
	{eqn='euler', gridSize={256,256}, initState='Sedov', solver='hll', integrator='forward Euler', fluxLimiter='superbee'},
	{eqn='euler', gridSize={256,256}, initState='Sedov', solver='hll', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
	{eqn='euler', gridSize={256,256}, initState='Sedov', solver='hll', integrator='Runge-Kutta 4, TVD', usePLM='plm-eig-prim-ref'},
}

for _,cfg in ipairs(configurations) do

	local args = {
		app = self,
		eqn = cfg.eqn,
		dim = #cfg.gridSize,
		integrator = cfg.integrator,
		fluxLimiter = cfg.fluxLimiter,
		usePLM = cfg.usePLM,
		geometry = 'cartesian',
		mins = {-1,-1,-1},
		maxs = {1,1,1},
		gridSize = cfg.gridSize,
		boundary = {
			xmin = 'freeflow',
			xmax = 'freeflow',
			ymin = 'freeflow',
			ymax = 'freeflow',
			zmin = 'freeflow',
			zmax = 'freeflow',
		},
		initState = cfg.initState,
	}
		
	local destMovieName = table{
		'eqn='..args.eqn,
		'solver='..cfg.solver,
		'integrator='..args.integrator,
	}:append(args.usePLM 
		-- plm:
		and table{
			'plm='..args.usePLM,
		}:append(
			args.slopeLimiter and {'slopeLimiter='..args.slopeLimiter} or nil
		) 
		-- non-plm: use flux limiter
		or (
			args.fluxLimiter and {'fluxLimiter='..args.fluxLimiter} or nil
		)
	):append{	
		'init='..args.initState,
		'gridSize='..range(args.dim):map(function(i) return args.gridSize[i] end):concat'x',
	}:concat', '..'.mp4'
	print(destMovieName)

	if io.fileexists(rundir..'/'..destMovieName) then
		print("I already found movie "..destMovieName)
	else

		local App = class(require 'app')
		function App:setup(clArgs)
			args.app = self
			self.solvers:insert(require('solver.'..cfg.solver)(args))
			sdl.SDL_SetWindowSize(self.window, 1280, 720)
			self.running = true
			self.exitTime = 1
			self.createAnimation = true
			self.displayVectorField_step = 8
		end
		local app  = App()
		
		app:run()

		-- once we're done ...
		local ext = app.screenshotExts[app.screenshotExtIndex]
		local dir = 'screenshots/'..app.screenshotDir
		assert(run('ffmpeg -y -i "'..dir..'/%05d.'..ext..'" "'..rundir..'/'..destMovieName..'"'))
		for i=0,app.screenshotIndex-1 do
			os.remove(dir..('/%05d.'):format(i)..ext)
		end
		local sep = ffi.os == 'Windows' and '\\' or '/'
		run('rmdir '..dir:gsub('/', sep))
	end
end

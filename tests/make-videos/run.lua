#!/usr/bin/env luajit
--[[
this will make movies of predefined configurations
it expects ffmpeg to be installed
--]]

local ffi = require 'ffi'
local unistd = require 'ffi.c.unistd'

require 'ffi.c.stdlib'
local rundirp = unistd.getcwd(nil, 0)
local rundir = ffi.string(rundirp)
ffi.C.free(rundirp)
unistd.chdir'../..'

local sdl = require 'ffi.sdl'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local os = require 'ext.os'
local file = require 'ext.file'

cmdline = {sys='console'}	-- set this before require 'hydro.app'

for k,v in pairs(require 'tests.util') do _G[k] = v end

local configurations = 
outer(
	{ {gridSize={256}} },
	-- [[ Euler & SRHD
	outer(
		{ {eqn='euler'} },
		{
			{initState='Sod', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
			{initState='Sod', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
			{initState='Sod', solver='roe', integrator='backward Euler', fluxLimiter='superbee'},
			{initState='Sod', solver='roe', integrator='backward Euler', usePLM='plm-eig-prim-ref'},
			
			{initState='Sod', solver='hll', integrator='forward Euler', fluxLimiter='superbee'},
			{initState='Sod', solver='hll', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
			{initState='Sod', solver='hll', integrator='backward Euler', fluxLimiter='superbee'},
			{initState='Sod', solver='hll', integrator='backward Euler', usePLM='plm-eig-prim-ref'},
			
			{initState='Sod', solver='euler-burgers', integrator='forward Euler', fluxLimiter='superbee'},
			{eqn='srhd', initState='Sod', solver='srhd-roe', integrator='forward Euler', fluxLimiter='superbee'},
			--{initState='Sod', solver='srhd-roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},		-- not working yet

			{initState='Sod', eqn='mhd', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
			{initState='Sod', eqn='mhd', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
			{initState='Sod', eqn='mhd', solver='roe', integrator='backward Euler', fluxLimiter='superbee'},
			{initState='Sod', eqn='mhd', solver='roe', integrator='backward Euler', usePLM='plm-eig-prim-ref'},


			{initState='Sedov', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
			{initState='Sedov', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
			{initState='Sedov', solver='roe', integrator='backward Euler', usePLM='plm-eig-prim-ref'},		-- b.e. + PLM gets nans
			{initState='Sedov', solver='hll', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
			{eqn='srhd', initState='Sedov', solver='srhd-roe', integrator='forward Euler', fluxLimiter='superbee'},
			--{eqn='srhd', initState='Sedov', solver='srhd-roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
			--{eqn='srhd', initState='Sedov', solver='srhd-roe', integrator='backward Euler', usePLM='plm-eig-prim-ref'},
	

			{initState='self-gravitation test 1', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
			{initState='self-gravitation test 1', solver='hll', integrator='forward Euler', fluxLimiter='superbee'},
			{eqn='srhd', initState='self-gravitation test 1', solver='srhd-roe', integrator='forward Euler', fluxLimiter='superbee'},
			--{eqn='srhd', initState='self-gravitation test 1', solver='hll', integrator='forward Euler', fluxLimiter='superbee'},
		}
	)
	--]]
	-- [[ SRHD
	:append(outer(
		{
			{eqn='srhd', solver='srhd-roe', integrator='forward Euler', fluxLimiter='superbee'},
			{eqn='srhd', solver='srhd-roe', integrator='backward Euler', fluxLimiter='superbee'},
		},
		{	
			--{initState='relativistic shock reflection'},			-- not working.  these initial conditions are constant =P
			{initState='relativistic blast wave test problem 1'},
			{initState='relativistic blast wave test problem 2'},
			{initState='relativistic blast wave interaction'},
		}
	))
	--]]
	-- [[ ideal MHD
	:append{
		{initState='Brio-Wu', eqn='mhd', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
		{initState='Brio-Wu', eqn='mhd', solver='roe', integrator='backward Euler', fluxLimiter='superbee'},
		{initState='Brio-Wu', eqn='mhd', solver='roe', integrator='forward Euler', usePLM='plm-eig'},
		{initState='Brio-Wu', eqn='mhd', solver='roe', integrator='backward Euler', usePLM='plm-eig'},
		
		--plm-eig-prim-ref not working with mhd...
		--{initState='Brio-Wu', eqn='mhd', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
		--{initState='Brio-Wu', eqn='mhd', solver='roe', integrator='backward Euler', usePLM='plm-eig-prim-ref'},
		-- also not working	
		--{initState='Brio-Wu', eqn='mhd', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim'},
		--{initState='Brio-Wu', eqn='mhd', solver='roe', integrator='backward Euler', usePLM='plm-eig-prim'},
		
	}
	--]]
	-- [[ Maxwell
	:append(outer(
		{
			{eqn='maxwell', solver='roe', integrator='forward Euler', fluxLimiter='superbee'}, 
			{eqn='maxwell', solver='roe', integrator='backward Euler', fluxLimiter='superbee'}, 
			{eqn='maxwell', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'}, 
			{eqn='maxwell', solver='roe', integrator='backward Euler', usePLM='plm-eig-prim-ref'}, 
			{eqn='maxwell', solver='hll', integrator='forward Euler', fluxLimiter='superbee'}, 
			{eqn='maxwell', solver='hll', integrator='backward Euler', fluxLimiter='superbee'}, 
			{eqn='maxwell', solver='hll', integrator='forward Euler', usePLM='plm-eig-prim-ref'}, 
			{eqn='maxwell', solver='hll', integrator='backward Euler', usePLM='plm-eig-prim-ref'}, 
		},
		{
			{initState='Maxwell default'},
			{initState='Maxwell scattering around cylinder', movieFrameDT=.01, movieEndTime=5},
			--{initState='Maxwell scattering around Koch snowflake'},	-- 2D only
			--{initState='Maxwell wire'},	-- 2D / 3D only
		}
	))
	--]]
	--[[ GR
	:append(outer(
		{
			{eqn='adm1d_v1', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
			{eqn='adm1d_v2', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
			{eqn='adm3d', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
			{eqn='bssnok-fd', solver='bssnok-fd', integrator='backward Euler'},
		},
		{
			{initState='gaussian perturbation'},
			{initState='plane gauge wave'},
			{initState='plane gauge wave'},
		}
	))
	--]]
):append(
	-- [[ Euler & SRHD
	outer(
	{
		{eqn='euler', gridSize={256,256} },
	},
	table{
		{initState='Sod', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
		{initState='Sod', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
		{initState='Sod', solver='roe', integrator='backward Euler', fluxLimiter='superbee'},
		{initState='Sod', solver='hll', integrator='forward Euler', fluxLimiter='superbee'},
		{initState='Sod', solver='hll', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
		{eqn='srhd', initState='Sod', solver='srhd-roe', integrator='forward Euler', fluxLimiter='superbee'},
		--{eqn='srhd', initState='Sod', solver='srhd-roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
		
		{initState='Sedov', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
		{initState='Sedov', solver='roe', integrator='backward Euler', fluxLimiter='superbee'},
		{initState='Sedov', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
		{initState='Sedov', solver='roe', integrator='Runge-Kutta 4, TVD', usePLM='plm-eig-prim-ref'},
		
		{initState='Sedov', solver='hll', integrator='forward Euler', fluxLimiter='superbee'},
		{initState='Sedov', solver='hll', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
		{initState='Sedov', solver='hll', integrator='Runge-Kutta 4, TVD', usePLM='plm-eig-prim-ref'},
	}:append(outer(
		{
			{solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
			{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
			{solver='roe', integrator='backward Euler', usePLM='plm-eig-prim-ref'},
			{solver='hll', integrator='forward Euler', fluxLimiter='superbee'},
			{solver='hll', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
			{solver='hll', integrator='backward Euler', usePLM='plm-eig-prim-ref'},
		},
		{
			{initState='Kelvin-Helmholtz', movieFrameDT=.1, movieEndTime=10},
			{initState='shock bubble interaction', movieEndTime=5},
			{initState='configuration 1', movieEndTime=.3},
			{initState='configuration 2', movieEndTime=.3},
			{initState='configuration 3', movieEndTime=.3},
			{initState='configuration 4', movieEndTime=.3},
			{initState='configuration 5', movieEndTime=.3},
			{initState='configuration 6', movieEndTime=.3},
		}
	)):append(outer(
		{
			{solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
			{solver='roe', integrator='backward Euler', fluxLimiter='superbee'},
			{solver='hll', integrator='forward Euler', fluxLimiter='superbee'},
			{eqn='srhd', solver='srhd-roe', integrator='forward Euler', fluxLimiter='superbee'},
			{eqn='srhd', solver='srhd-roe', integrator='backward Euler', fluxLimiter='superbee'},
		},
		{
			{initState='self-gravitation test 1'},
			{initState='self-gravitation test 1 spinning'},
			{initState='self-gravitation test 2'},
			{initState='self-gravitation test 2 orbiting'},
		}
	))
	--]]
	-- [[ SRHD
	:append(
		outer(
			{
				{eqn='srhd', solver='srhd-roe', integrator='forward Euler', fluxLimiter='superbee'},
				{eqn='srhd', solver='srhd-roe', integrator='backward Euler', fluxLimiter='superbee'},
			},
			{	
				{initState='relativistic shock reflection'},			-- not working.  these initial conditions are constant =P
				{initState='relativistic blast wave test problem 1'},
				{initState='relativistic blast wave test problem 2'},
				{initState='relativistic blast wave interaction'},
			}
		)
	)
	--]]
	-- [[ ideal MHD
	:append{
		{initState='Orszag-Tang', eqn='mhd', solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
		{initState='Orszag-Tang', eqn='mhd', solver='roe', integrator='backward Euler', fluxLimiter='superbee'},
		{initState='Orszag-Tang', eqn='mhd', solver='roe', integrator='forward Euler', usePLM='plm-eig'},
		{initState='Orszag-Tang', eqn='mhd', solver='roe', integrator='backward Euler', usePLM='plm-eig'},
	}
	--]]
	-- [[ Maxwell
	:append(outer(
		{
			{eqn='maxwell', solver='roe', integrator='forward Euler', fluxLimiter='superbee'}, 
			{eqn='maxwell', solver='roe', integrator='backward Euler', fluxLimiter='superbee'}, 
			{eqn='maxwell', solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'}, 
			{eqn='maxwell', solver='roe', integrator='backward Euler', usePLM='plm-eig-prim-ref'}, 
			{eqn='maxwell', solver='hll', integrator='forward Euler', fluxLimiter='superbee'}, 
			{eqn='maxwell', solver='hll', integrator='backward Euler', fluxLimiter='superbee'}, 
			{eqn='maxwell', solver='hll', integrator='forward Euler', usePLM='plm-eig-prim-ref'}, 
			{eqn='maxwell', solver='hll', integrator='backward Euler', usePLM='plm-eig-prim-ref'}, 
		},
		{
			{initState='Maxwell default'},
			{initState='Maxwell scattering around cylinder', movieFrameDT=.01, movieEndTime=5},
			{initState='Maxwell scattering around Koch snowflake', movieFrameDT=.01, movieEndTime=5},	-- 2D only
			{initState='Maxwell wire'},	-- 2D / 3D only
		}
	))
	--]]
))

for _,cfg in ipairs(configurations) do

	local args = {
		app = self,
		eqn = cfg.eqn,
		dim = #cfg.gridSize,
		integrator = cfg.integrator,
		fluxLimiter = cfg.fluxLimiter,
		usePLM = cfg.usePLM,
		coord = 'cartesian',
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
		
	local destMovieName = nameForConfig(cfg, args)..'.mp4'
	print(destMovieName)
	
	if file(rundir..'/'..destMovieName):exists() then
		print("I already found movie "..destMovieName)
	else
		local movieStartTime = cfg.movieStartTime or 0
		local movieEndTime = cfg.movieEndTime or 1
		local movieFrameDT = cfg.movieFrameDT or 0 

		local App = class(require 'hydro.app')
		function App:setup(clArgs)
			args.app = self
			self.solvers:insert(require('hydro.solver.'..cfg.solver)(args))
			sdl.SDL_SetWindowSize(self.window, 1280, 720)
			
			-- why am I getting weird graphics glitches upon startup?
			-- first it was imgui giving a DeltaTime>=0 assert fail only for the *second* app run in a batch
			-- now it is the apps losing all their textures ... which I think might be SDL_SetWindowSize 
			sdl.SDL_Delay(1)
			
			self.running = true
			self.exitTime = movieEndTime + movieStartTime
			self.displayVectorField_step = 8
		end
		
		local recording
		local nextCaptureTime
		function App:update(...)
	
			local oldestSolver = self.solvers:inf(function(a,b) return a.t < b.t end)
			
			-- if the recording hasn't started then start it ... 
			-- ... as soon as it passes the movieStartTime 
			if not recording then
				if not movieStartTime 
				or oldestSolver.t >= movieStartTime
				then
					recording = true
				end
			end

			if recording then
				if not nextCaptureTime
				or nextCaptureTime <= oldestSolver.t
				then
					self.createAnimation = 'once'
					nextCaptureTime = oldestSolver.t + movieFrameDT
				end
			end

			return App.super.update(self, ...)
		end
		local app  = App()
		
		app:run()

		-- once we're done ...
		local ext = app.screenshotExts[app.screenshotExtIndex]
		local dir = 'screenshots/'..app.screenshotDir
		assert(run('ffmpeg -y -i "'..dir..'/%05d.'..ext..'" "'..rundir..'/'..destMovieName..'"'))
		for i=0,app.screenshotIndex-1 do
			file(dir..('/%05d.'):format(i)..ext):remove()
		end
		run('rmdir '..dir:gsub('/', os.sep))
--for some odd reason, every time it runs, it messes up the GL state ...
os.exit(1)
	end
end

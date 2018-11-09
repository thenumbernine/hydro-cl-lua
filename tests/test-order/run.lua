#!/usr/bin/env luajit
--[[
this will run some order of accuracy tests on different configurations 
--]]

local ffi = require 'ffi'

local clnumber = require 'cl.obj.number'
local template = require 'template'
local sdl = require 'ffi.sdl'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local fromlua = require 'ext.fromlua'
local tolua = require 'ext.tolua'
local file = require 'ext.file'
local string = require 'ext.string'
local io = require 'ext.io'
local matrix = require 'matrix'
local gnuplot = require 'gnuplot'

local rundir = string.trim(io.readproc'pwd')

-- from here on the require's expect us to be in the hydro-cl directory
-- I should change this, and prefix all hydro-cl's require()s with 'hydro-cl', so it is require()able from other projects
require 'ffi.c.unistd'
ffi.C.chdir'../..'

__useConsole__ = true	-- set this before require 'app'

for k,v in pairs(require 'tests.util') do _G[k] = v end

local configurations = outer(
	{
		{
			eqn='euler',
			initState = 'advect wave',
		}
	},
		-- final error at n=1024 on the right:
	{	
		-- schemes
		--{solver='weno5', integrator='forward Euler'},									-- 0.099267582810394
		--{solver='hll', integrator='forward Euler'},									-- 0.00060136599076404
		--{solver='euler-hllc', integrator='forward Euler'},							-- 0.00048873499978618
		--{solver='euler-burgers', integrator='forward Euler'},							-- 0.0004752949543945

		-- flux-limiters w/roe scheme:
		--{solver='roe', integrator='forward Euler', fluxLimiter='smart'},				-- fails on n=1024
		--{solver='roe', integrator='forward Euler', fluxLimiter='ospre'},				-- 0.99999990000004
		--{solver='roe', integrator='forward Euler', fluxLimiter='Fromm'},				-- 0.99999990000003
		--{solver='roe', integrator='forward Euler', fluxLimiter='CHARM'},				-- 0.99999990000003
		--{solver='roe', integrator='forward Euler', fluxLimiter='van Albada 1'},		-- 0.99999990000003
		--{solver='roe', integrator='forward Euler', fluxLimiter='Barth-Jespersen'},	-- 0.99999990000003
		--{solver='roe', integrator='forward Euler', fluxLimiter='Beam-Warming'},		-- 0.99999990000002
		--{solver='roe', integrator='forward Euler', fluxLimiter='van Leer'},			-- 0.99999990000017
		--{solver='roe', integrator='forward Euler', fluxLimiter='van Albada 2'},		-- 0.9999999
		--{solver='roe', integrator='forward Euler', fluxLimiter='donor cell'},			-- 0.00048873499978677
		--{solver='roe', integrator='forward Euler', fluxLimiter='Oshker'},				-- 8.1594383698357e-05
		--{solver='roe', integrator='forward Euler', fluxLimiter='superbee'},			-- 8.0220379351244e-05
		--{solver='roe', integrator='forward Euler', fluxLimiter='Sweby'},				-- 7.5406302510808e-05
		--{solver='roe', integrator='forward Euler', fluxLimiter='HQUICK'},				-- 7.3985100798005e-05
		--{solver='roe', integrator='forward Euler', fluxLimiter='Koren'},				-- 7.3374191905257e-05
		--{solver='roe', integrator='forward Euler', fluxLimiter='HCUS'},				-- 7.155162562564e-05
		--{solver='roe', integrator='forward Euler', fluxLimiter='monotized central'},	-- 6.8474331334665e-05
		--{solver='roe', integrator='forward Euler', fluxLimiter='UMIST'},				-- 6.3705455038239e-05
		--{solver='roe', integrator='forward Euler', fluxLimiter='minmod'},				-- 5.9129797191892e-05
		{solver='roe', integrator='forward Euler', fluxLimiter='Lax-Wendroff'},			-- 1.2048891136515e-06

		-- my PLM attempts:
		--{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},			-- 0.00049148119638364
		--{solver='roe', integrator='forward Euler', usePLM='plm-eig'},						-- 0.00049103769947135
		--{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim'},				-- 0.00049075156058252
		--{solver='roe', integrator='forward Euler', usePLM='plm-cons'},					-- 0.00048873499978677
		--{solver='roe', integrator='forward Euler', usePLM='plm-athena'},					-- 0.00014403561237557

		-- various explicit integrators:
		--{solver='roe', integrator='Runge-Kutta 2, non-TVD', fluxLimiter='Lax-Wendroff'},	-- 0.12722668099294
		--{solver='roe', integrator='Runge-Kutta 2 Heun', fluxLimiter='Lax-Wendroff'},		-- 0.090030643072756
		--{solver='roe', integrator='Runge-Kutta 2, TVD', fluxLimiter='Lax-Wendroff'},		-- 0.090030643072756
		--{solver='roe', integrator='Runge-Kutta 3', fluxLimiter='Lax-Wendroff'},			-- 0.08997529488926
		--{solver='roe', integrator='Runge-Kutta 3, TVD', fluxLimiter='Lax-Wendroff'},		-- 0.063645931031262
		--{solver='roe', integrator='Runge-Kutta 2 Ralston', fluxLimiter='Lax-Wendroff'},	-- 0.048729889553282
		--{solver='roe', integrator='Runge-Kutta 4', fluxLimiter='Lax-Wendroff'},			-- 0.032939577524704
		--{solver='roe', integrator='Runge-Kutta 4, TVD', fluxLimiter='Lax-Wendroff'},		-- 0.032939577524698
		--{solver='roe', integrator='Runge-Kutta 4, non-TVD', fluxLimiter='Lax-Wendroff'},	-- 0.03293957752469
		--{solver='roe', integrator='Runge-Kutta 4, 3/8ths rule', fluxLimiter='Lax-Wendroff'},-- 0.024844936602353
		--{solver='roe', integrator='Runge-Kutta 2', fluxLimiter='Lax-Wendroff'},			-- 1.2048891136515e-06

		-- implicit integrators:
		-- backward euler with epsilon=1e-10
		--{solver='weno5', integrator='backward Euler'},									-- 0.09925195780133
		
		-- all of these dip down at some optimal size for their epsilon, then pop back up
		--{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-10}},	-- 9.4498933891175e-06
		--{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-20}},	-- 7.0082322822353e-06
		-- except this one, which just follows forward-Euler exactly
		--{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-30}},	-- 1.178306878477e-06
		
		{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-10}},	-- 9.4498933891175e-06
		{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-20}},	-- 8.1079737186957e-06
		{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-30}},	-- 1.1257993855623e-06
	}
)

local dim = 1
local sizes = range(3,10):map(function(x) return 2^x end)
local errorsForConfig = table()
local errorNames = table()
for _,cfg in ipairs(configurations) do

	local args = {
		app = self,
		cfl = .5,
		eqn = cfg.eqn,
		dim = dim,
		integrator = cfg.integrator,
		integratorArgs = cfg.integratorArgs,
		fluxLimiter = cfg.fluxLimiter,
		usePLM = cfg.usePLM,
		coord = 'cartesian',
		mins = {-1,-1,-1},
		maxs = {1,1,1},
		initState = cfg.initState,
	}
	
	local destName = table{
		'solver='..cfg.solver,
		'integrator='..args.integrator:gsub('/','-'),
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
	):append(
		args.integrator == 'backward Euler' 
		and args.integratorArgs 
		and table{
			args.integratorArgs.epsilon and ('be.eps='..args.integratorArgs.epsilon) or nil
		}:append{
			args.integratorArgs.restart and ('be.restart='..args.integratorArgs.restart) or nil
		} or nil
	):concat', '
	print(destName)

	local errors 
	local srcfn = rundir..'/'..destName..'.lua'
	local srcdata = file[srcfn]
	if srcdata then
		local res = fromlua(srcdata)
		errors = res.errors
		local storedSizes = res.sizes
		if matrix(sizes) ~= matrix(storedSizes) then
			errors = nil
		end
	end
	if not errors then
		errors = table()
		for _,size in ipairs(sizes) do
			args.gridSize = {size}
	print()
	print(size)
	print()
			
			-- TODO make sure solver_t->init_v0x == 1/duration and solver_t->maxs.x - mins.x == 2
			-- otherwise, for durations t=100 and t=1 the results look close enough to the same
			local duration = cfg.duration or 1
			
			local App = class(require 'app')
			function App:setup(clArgs)
				args.app = self
				local solver = require('solver.'..cfg.solver)(args)
				self.solvers:insert(solver)
				self.exitTime = duration
				self.running = true
			end
			function App:requestExit()
				App.super.requestExit(self)
			
				-- now compare the U buffer to the exact 
				assert(#self.solvers == 1)
				local solver = self.solvers[1]
				local _, var = solver.displayVars:find(nil, function(var) return var.name == 'U rho' end)
				assert(var, "failed to find U rho var")
				solver:calcDisplayVarToBuffer(var)	
				-- now in solver.reduceBuf
				local ptr = solver.calcDisplayVarToTexPtr
				local numCells = solver.numCells
				local numGhost = solver.numGhost
				self.cmds:enqueueReadBuffer{buffer=solver.reduceBuf, block=true, size=ffi.sizeof(self.real) * numCells, ptr=ptr}
				-- now in ptr
				
				local xmin = solver.solverPtr.mins.x
				local xmax = solver.solverPtr.maxs.x
				local width = xmax - xmin
				local xs = range(numCells-2*numGhost):mapi(function(i)
					return (i-.5) * solver.solverPtr.grid_dx.x + solver.solverPtr.mins.x
				end)
				local ys = range(numCells-2*numGhost):mapi(function(i)
					return ptr[i+numGhost-1]
				end)
				local rho0 = solver.solverPtr.init_rho0
				local rho1 = solver.solverPtr.init_rho1
				local u0 = solver.solverPtr.init_v0x
				local k0 = 2 * math.pi / width
				local t1 = solver.t
				local exact = xs:map(function(x,i)
					return rho0 + rho1 * math.cos(k0 * (x - u0 * t1))
				end)

				local diff = matrix(exact) - matrix(ys)
				errors:insert(diff:normL1() / #diff)
			
				--[[ plotting immediately
				gnuplot{
					output = rundir..'/'..destName,
					style = 'data lines',
					data = {xs, ys, exact},
					{using = '1:2', title=''..size},
					{using = '1:3', title='exact'},
				}
				--]]
			end
			
			local app  = App()
			app:run()
		
			if file.stop then file.stop = nil os.exit(1) end
		end
		file[srcfn] = tolua{errors=errors, sizes=sizes}
	end
	
	errorsForConfig:insert(errors)
	errorNames:insert(destName)
end

-- [[ plot errors 
gnuplot(
	table({
		output = rundir..'/results.png',
		terminal = 'png size 1200,700',
		style = 'data linespoints',
		log = 'xy',
		xlabel = 'grid size',
		ylabel = 'L1 error',
		data = table{sizes}:append(errorsForConfig),
	},
	errorNames:mapi(function(name,i)
		return {using = '1:'..(i+1), title=name}
	end)
))
--]]


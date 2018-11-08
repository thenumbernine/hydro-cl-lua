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
			initState = 'advect wave',
		},
	},
	outer(
		{ {eqn='euler'} },
		{
			{solver='euler-burgers', integrator='forward Euler'},
--			{solver='hll', integrator='forward Euler'},
--			{solver='euler-hllc', integrator='forward Euler'},
			--{solver='roe', integrator='forward Euler', fluxLimiter='donor cell'},
			--{solver='roe', integrator='forward Euler', fluxLimiter='Lax-Wendroff'},
			--{solver='roe', integrator='forward Euler', fluxLimiter='minmod'},
			--{solver='roe', integrator='forward Euler', fluxLimiter='monotized central'},
--			{solver='roe', integrator='forward Euler', fluxLimiter='superbee'},
			
--			{solver='weno5', integrator='forward Euler'}, -- usePLM='plm-cons'
			
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons'},
			--{solver='roe', integrator='forward Euler', usePLM='plm-eig'},
			--{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim'},
--			{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},
			--{solver='roe', integrator='forward Euler', usePLM='plm-athena'},
			
			--{solver='roe', usePLM='plm-eig-prim-cons', integrator = 'Runge-Kutta 2'},
			--{solver='roe', usePLM='plm-eig-prim-cons', integrator = 'Runge-Kutta 2 Heun'},
			--{solver='roe', usePLM='plm-eig-prim-cons', integrator = 'Runge-Kutta 2 Ralston'},
			--{solver='roe', usePLM='plm-eig-prim-cons', integrator = 'Runge-Kutta 3'},
			--{solver='roe', usePLM='plm-eig-prim-cons', integrator = 'Runge-Kutta 4'},
			--{solver='roe', usePLM='plm-eig-prim-cons', integrator = 'Runge-Kutta 4, 3/8ths rule'},
			--{solver='roe', usePLM='plm-eig-prim-cons', integrator = 'Runge-Kutta 2, TVD'},
			--{solver='roe', usePLM='plm-eig-prim-cons', integrator = 'Runge-Kutta 2, non-TVD'},
			--{solver='roe', usePLM='plm-eig-prim-cons', integrator = 'Runge-Kutta 3, TVD'},
--			{solver='roe', usePLM='plm-eig-prim-cons', integrator = 'Runge-Kutta 4, TVD'},
			--{solver='roe', usePLM='plm-eig-prim-cons', integrator = 'Runge-Kutta 4, non-TVD'},

--			{solver='roe', integrator='backward Euler', usePLM='plm-eig-prim-ref'},
--			{solver='weno5', integrator='backward Euler'}, --usePLM='plm-eig-prim-ref'
		}
	)
)

local dim = 1
local sizes = range(3,10):map(function(x) return 2^x end)
local errorsForConfig = table()
local errorNames = table()
for _,cfg in ipairs(configurations) do

	local args = {
		app = self,
		eqn = cfg.eqn,
		dim = dim,
		integrator = cfg.integrator,
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
	):concat', '
	print(destName)

	local errors 
	local srcfn = rundir..'/'..destName..'.lua'
	local srcdata = file[srcfn]
	if srcdata then
		errors = fromlua(srcdata)
		if matrix(table.keys(errors):sort()) ~= matrix(sizes) then
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
			local startTime = cfg.startTime or 0
			local duration = cfg.duration or 1
			
			local App = class(require 'app')
			function App:setup(clArgs)
				args.app = self
				local solver = require('solver.'..cfg.solver)(args)
				self.solvers:insert(solver)
				self.exitTime = duration + startTime
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
					return rho0 + rho1 * math.sin(k0 * (x - u0 * t1))
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
		file[srcfn] = tolua(errors)
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


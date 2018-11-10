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
local math = require 'ext.math'
local file = require 'ext.file'
local string = require 'ext.string'
local io = require 'ext.io'
local matrix = require 'matrix'
local gnuplot = require 'gnuplot'
require 'ffi.c.unistd'

local rundir = string.trim(io.readproc'pwd')

-- from here on the require's expect us to be in the hydro-cl directory
-- I should change this, and prefix all hydro-cl's require()s with 'hydro-cl', so it is require()able from other projects
ffi.C.chdir'../..'

for k,v in pairs(require 'tests.util') do _G[k] = v end

__useConsole__ = true	-- set this before require 'app'


local problems = {}

problems['advect wave'] = {
	configurations = outer(
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
			{solver='hll', integrator='forward Euler'},									-- 0.00060136599076404
			{solver='euler-hllc', integrator='forward Euler'},							-- 0.00048873499978618
			{solver='euler-burgers', integrator='forward Euler'},							-- 0.0004752949543945

			-- flux-limiters w/roe scheme:
			{solver='roe', integrator='forward Euler', fluxLimiter='smart'},				-- fails on n=1024
			{solver='roe', integrator='forward Euler', fluxLimiter='ospre'},				-- 0.99999990000004
			{solver='roe', integrator='forward Euler', fluxLimiter='Fromm'},				-- 0.99999990000003
			{solver='roe', integrator='forward Euler', fluxLimiter='CHARM'},				-- 0.99999990000003
			{solver='roe', integrator='forward Euler', fluxLimiter='van Albada 1'},		-- 0.99999990000003
			{solver='roe', integrator='forward Euler', fluxLimiter='Barth-Jespersen'},	-- 0.99999990000003
			{solver='roe', integrator='forward Euler', fluxLimiter='Beam-Warming'},		-- 0.99999990000002
			{solver='roe', integrator='forward Euler', fluxLimiter='van Leer'},			-- 0.99999990000017
			{solver='roe', integrator='forward Euler', fluxLimiter='van Albada 2'},		-- 0.9999999
			{solver='roe', integrator='forward Euler', fluxLimiter='donor cell'},			-- 0.00048873499978677
			{solver='roe', integrator='forward Euler', fluxLimiter='Oshker'},				-- 8.1594383698357e-05
			{solver='roe', integrator='forward Euler', fluxLimiter='superbee'},			-- 8.0220379351244e-05
			{solver='roe', integrator='forward Euler', fluxLimiter='Sweby'},				-- 7.5406302510808e-05
			{solver='roe', integrator='forward Euler', fluxLimiter='HQUICK'},				-- 7.3985100798005e-05
			{solver='roe', integrator='forward Euler', fluxLimiter='Koren'},				-- 7.3374191905257e-05
			{solver='roe', integrator='forward Euler', fluxLimiter='HCUS'},				-- 7.155162562564e-05
			{solver='roe', integrator='forward Euler', fluxLimiter='monotized central'},	-- 6.8474331334665e-05
			{solver='roe', integrator='forward Euler', fluxLimiter='UMIST'},				-- 6.3705455038239e-05
			{solver='roe', integrator='forward Euler', fluxLimiter='minmod'},				-- 5.9129797191892e-05
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
			
			--{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-10}},	-- 9.4498933891175e-06
			--{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-20}},	-- 8.1079737186957e-06
			--{solver='roe', fluxLimiter='Lax-Wendroff', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-30}},	-- 1.1257993855623e-06
		}
	),
	makeExact = function(problem, xs, solver)
		local rho0 = solver.solverPtr.init_rho0
		local rho1 = solver.solverPtr.init_rho1
		local u0 = solver.solverPtr.init_v0x
		local xmin = solver.solverPtr.mins.x
		local xmax = solver.solverPtr.maxs.x
		local width = xmax - xmin
		local k0 = 2 * math.pi / width
		local t1 = solver.t
		return xs:map(function(x,i)
			return rho0 + rho1 * math.cos(k0 * (x - u0 * t1))
		end)
	end,

	-- TODO make sure solver_t->init_v0x == 1/duration and solver_t->maxs.x - mins.x == 2
	-- otherwise, for durations t=100 and t=1 the results look close enough to the same
	duration = 1,
}


problems.Sod = {
	-- copy of the above problem ... maybe put somewhere else
	configurations = outer(
		{
			{
				eqn='euler',
				initState = 'Sod',
			}
		},
			-- final error at n=1024 on the right:
		{	
			-- schemes
			--{solver='weno5', integrator='forward Euler'},									-- 0.028050334485117
			--{solver='hll', integrator='forward Euler'},									-- 0.0039886633966807
			--{solver='euler-hllc', integrator='forward Euler'},							-- 0.0036984733332097
			--{solver='euler-burgers', integrator='forward Euler'},							-- 0.0031694650615551

			-- flux-limiters w/roe scheme:
			--{solver='roe', integrator='forward Euler', fluxLimiter='Lax-Wendroff'},		-- 0.47968895872093
			--{solver='roe', integrator='forward Euler', fluxLimiter='donor cell'},			-- 0.0036660453357688
			--{solver='roe', integrator='forward Euler', fluxLimiter='ospre'},				-- 0.0030672192288732
			--{solver='roe', integrator='forward Euler', fluxLimiter='Beam-Warming'},		-- 0.0022967397811603
			--{solver='roe', integrator='forward Euler', fluxLimiter='minmod'},				-- 0.0020333322754062
			--{solver='roe', integrator='forward Euler', fluxLimiter='van Albada 2'},		-- 0.0020307570851758
			--{solver='roe', integrator='forward Euler', fluxLimiter='UMIST'},				-- 0.0019615353080538
			--{solver='roe', integrator='forward Euler', fluxLimiter='Oshker'},				-- 0.0019441809134923
			--{solver='roe', integrator='forward Euler', fluxLimiter='CHARM'},				-- 0.0019170086674662
			--{solver='roe', integrator='forward Euler', fluxLimiter='van Leer'},			-- 0.0019121788850881
			--{solver='roe', integrator='forward Euler', fluxLimiter='HQUICK'},				-- 0.001908042213926
			--{solver='roe', integrator='forward Euler', fluxLimiter='Fromm'},				-- 0.0019069360974506
			--{solver='roe', integrator='forward Euler', fluxLimiter='Koren'},				-- 0.0018974481449319
			--{solver='roe', integrator='forward Euler', fluxLimiter='HCUS'},				-- 0.0018945018297761
			--{solver='roe', integrator='forward Euler', fluxLimiter='monotized central'},	-- 0.0018898857563959
			--{solver='roe', integrator='forward Euler', fluxLimiter='Barth-Jespersen'},	-- 0.0018828276652798
			{solver='roe', integrator='forward Euler', fluxLimiter='superbee'},			-- 0.00187534979448
			--{solver='roe', integrator='forward Euler', fluxLimiter='Sweby'},				-- 0.0018687046583677
			--{solver='roe', integrator='forward Euler', fluxLimiter='van Albada 1'},		-- 0.001865866089593
			--{solver='roe', integrator='forward Euler', fluxLimiter='smart'},				-- 0.0018130273610755	-- even though this is the lowest error of the flux limiters, it has a definite hiccup in it
			
			-- my PLM attempts:
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons'},				-- 0.0036660453357688
			--{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},		-- 0.0023632029528882	-- plm-eig-prim-ref is worse than plm-eig-prim
			--{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim'},			-- 0.0021682982996781	-- plm-eig-prim is worse than plm-eig
			--{solver='roe', integrator='forward Euler', usePLM='plm-eig'},					-- 0.0019239440708115	
			{solver='roe', integrator='forward Euler', usePLM='plm-athena'},				-- 0.0016947738688964	

			-- various explicit integrators:
			-- RK2 is exactly the same as forward Euler ... which means it's broken
			-- and everything else is worse
			--{solver='roe', integrator='Runge-Kutta 2, non-TVD', fluxLimiter='superbee'},		-- 0.30989158592001
			--{solver='roe', integrator='Runge-Kutta 3', fluxLimiter='superbee'},				-- 0.028049852882381
			--{solver='roe', integrator='Runge-Kutta 2 Heun', fluxLimiter='superbee'},			-- 0.0277197399486
			--{solver='roe', integrator='Runge-Kutta 2, TVD', fluxLimiter='superbee'},			-- 0.0277197399486
			--{solver='roe', integrator='Runge-Kutta 3, TVD', fluxLimiter='superbee'},			-- 0.019060751585123
			--{solver='roe', integrator='Runge-Kutta 2 Ralston', fluxLimiter='superbee'},		-- 0.014905909835352
			--{solver='roe', integrator='Runge-Kutta 4', fluxLimiter='superbee'},				-- 0.010605564924952
			--{solver='roe', integrator='Runge-Kutta 4, non-TVD', fluxLimiter='superbee'},		-- 0.010605564924952
			--{solver='roe', integrator='Runge-Kutta 4, TVD', fluxLimiter='superbee'},			-- 0.010605515773635
			--{solver='roe', integrator='Runge-Kutta 4, 3/8ths rule', fluxLimiter='superbee'},	-- 0.0083893332244032
			--{solver='roe', integrator='Runge-Kutta 2', fluxLimiter='superbee'},				-- 0.00187534979448	-- this is exactly the same as Forward Euler, so... something is wrong

			-- implicit integrators:
			-- backward euler with epsilon=1e-10
			--{solver='weno5', integrator='backward Euler'},									-- 0.028022560498779
		
			-- well, lower epsilon does better than higher epsilon 
			-- restart doesn't matter
			-- and implicit can't beat plm-athena
			--{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-30}},	-- 0.0018753132626329
			--{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-30}},	-- 0.0018703319992672
			--{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-20}},	-- 0.0018700680964215
			--{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-20}},	-- 0.0018586376742315
			--{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=10, epsilon=1e-10}},	-- 0.0017705698680532
			--{solver='roe', fluxLimiter='superbee', integrator='backward Euler', integratorArgs={restart=20, epsilon=1e-10}},	-- 0.0017705698680532
		
		}
	),

	-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
	makeExact = function(problem, xs, solver)
		-- TODO initial condition object to share these values with initialization
		local rhoL = solver.solverPtr.init_rhoL
		local rhoR = solver.solverPtr.init_rhoR
		local PL = solver.solverPtr.init_PL
		local PR = solver.solverPtr.init_PR
		local vL = 0
		local vR = 0
		local gamma = solver.solverPtr.heatCapacityRatio
		local t = solver.t

-- why
t = t * 4 / 3
		
		local muSq = (gamma - 1)/(gamma + 1)
		local K = PL / rhoL^gamma

		local CsL = math.sqrt(gamma * PL / rhoL)
		local CsR = math.sqrt(gamma * PR / rhoR)

		local solveP3 = function()
			local symmath = require 'symmath'
			local P3 = symmath.var'P3'
			local f = -2*CsL*(1 - (P3/PL)^((-1 + gamma)/(2*gamma)))/(CsR*(-1 + gamma)) + (-1 + P3/PR)*((1 - muSq)/(gamma*(muSq + P3/PR)))^.5
			local df_dP3 = f:diff(P3)()	
			local f_func = f:compile{P3}
			local df_dP3_func = df_dP3:compile{P3}
			local P3 = .5 * (PL + PR)
			local epsilon = 1e-16	-- this is the limit for the sod.f before it oscillates
			while true do
				local dP3 = -f_func(P3) / df_dP3_func(P3)
				--print(P3, dP3)
				if math.abs(dP3) <= epsilon then break end
				if not math.isfinite(dP3) then error('delta is not finite! '..tostring(dP3)) end
				P3 = P3 + dP3 
			end
			return P3
		end

		local P3 = solveP3()
		local P4 = P3

		local rho3 = rhoL * (P3 / PL) ^ (1 / gamma)

		local v3 = vR + 2 * CsL / (gamma - 1) * (1 - (P3 / PL)^((gamma - 1)/(2*gamma)))
		local v4 = v3

		local rho4 = rhoR * (P4 + muSq * PR) / (PR + muSq * P4)

		local vshock = v4 * rho4 / (rho4 - rhoR)
		local vtail = CsL - v4 / (1 - muSq)

		local v2 = function(x) return (1 - muSq) * (x/t + CsL) end
		-- Dullemon:
		--local rho2 = function(x) return (rhoL^gamma / (gamma * PL) * (v2(x) - x/t)^2)^(1/(gamma-1)) end
		-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
		local rho2 = function(x) return rhoL * (-muSq * (x / (CsL * t)) + (1 - muSq))^(2/(gamma-1)) end

		-- Dullemon:
		--local P2 = K * rho2^gamma
		-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
		local P2 = function(x) return PL * (-muSq * (x / (CsL * t)) + (1 - muSq)) ^ (2*gamma/(gamma-1)) end

		local C = function(c) return function() return c end end
		local function makefunc(t)
			return table(t):map(function(ti,i)
				return type(ti) == 'function'and t[i]or C(t[i])
			end)
		end

		local rhos = makefunc{rhoL, rho2, rho3, rho4, rhoR}
		local vs = makefunc{vL, v2, v3, v4, vR}
		local Ps = makefunc{PL, P2, P3, P4, PR}
		
		-- between regions 1 and 2
		local s1 = -CsL	

		-- between regions 2 and 3
		-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
		local s2 = -vtail

		local s3 = v3	-- between regions 3 and 4

		-- between regions 4 and 5 ...
		local s4 = vshock

		--print('wavespeeds:',s1,s2,s3,s4)

		local region = function(x)
			local xi = x / t
			if xi < s1 then
				return 1
			elseif xi < s2 then
				return 2
			elseif xi < s3 then
				return 3
			elseif xi < s4 then
				return 4
			else
				return 5
			end
		end

		return xs:map(function(x)
			local r = region(x)
			return rhos[r](x)	-- vs[r], Ps[r]
		end)
	end,
	
	duration = .2,
}

local problem = problems.Sod
--local problem = problems['advect wave']


local args = table{...}
local nocache = args:find'nocache'
local plotCompare = args:find'compare'


local dim = 1
local sizes = plotCompare and {1024} or range(3,10):map(function(x) return 2^x end)
local errorsForConfig = table()
local errorNames = table()
for _,cfg in ipairs(problem.configurations) do

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
	):append{
		'init='..args.initState,
	}:concat', '
	print(destName)

	local errors 
	local srcfn = rundir..'/'..destName..'.lua'
	local srcdata
	if not nocache then 
		srcdata = file[srcfn]
	end
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
			
			local duration = tonumber(problem.duration) or error("expected problem.duration")
			
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
				local exact = problem:makeExact(xs, solver)

				local diff = matrix(exact) - matrix(ys)
				errors:insert(diff:normL1() / #diff)
			
				if plotCompare then -- plotting immediately
					gnuplot{
						output = rundir..'/compare-graphs.png',
						style = 'data lines',
						data = {xs, ys, exact},
						{using = '1:2', title=''..size},
						{using = '1:3', title='exact'},
					}
					os.exit()				
				end
			end
			
			local app  = App()
			app:run()
		
			if file.stop then file.stop = nil os.exit(1) end
		end
		file[srcfn] = tolua{errors=errors, sizes=sizes}
	end
print(table.last(errors))	
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
		data = table{
			sizes:map(function(x) return 2/x end),	-- this assumes the domain is -1,1
		}:append(errorsForConfig),
	},
	errorNames:mapi(function(name,i)
		return {using = '1:'..(i+1), title=name}
	end)
))
--]]


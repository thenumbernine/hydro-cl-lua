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


local cmdline = {}
for _,w in ipairs(arg or {}) do
	local k,v = w:match'^(.-)=(.*)$'
	if k then
		cmdline[k] = fromlua(v)
	else
		cmdline[w] = true
	end
end

-- which problem to use
local problemName = cmdline.init or 'advect wave'
--local problemName = cmdline.init or 'Sod'

-- don't use cached results <-> regenerate results for selected tests
local nocache = cmdline.nocache

-- for the first configuration, run it at highest resolution, plot it with the exact solution, and quit
local plotCompare = cmdline.compare

-- for the first configuration, plot error vs time from 0 to duration
local plotErrorHistory = cmdline.history

-- exclusive with 'compare': don't use exact, instead use exponential regression
local uselin = cmdline.uselin


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
			--{solver='weno5', weno5method='1996 Jiang Shu', integrator='forward Euler'},		-- 0.092113212858329
			--{solver='weno5', weno5method='2008 Borges', integrator='forward Euler'},		-- 0.16680010146205	
			{solver='weno5', weno5method='2010 Shen Zha', integrator='forward Euler'},		-- 0.12853659670964	

			--{solver='hll', integrator='forward Euler'},									-- 0.00060136599076404
			--{solver='euler-burgers', integrator='forward Euler'},							-- 0.0004752949543945	

			-- why is RK3-TVD worse than forward Euler in all my hllc solvers?
			-- hllcMethod == 0
			--{solver='euler-hllc', hllcMethod=0, integrator='forward Euler', usePLM='plm-cons-alone'},   	-- 0.00060495293676367
			--{solver='euler-hllc', hllcMethod=0, integrator='forward Euler', usePLM='plm-prim'},         	-- 0.00060393770462868
			--{solver='euler-hllc', hllcMethod=0, integrator='forward Euler'},								-- 0.00048873499978739
			--{solver='euler-hllc', hllcMethod=0, integrator='forward Euler', usePLM='plm-cons'},         	-- 0.00048873499978739
			--{solver='euler-hllc', hllcMethod=0, integrator='Runge-Kutta 3, TVD', usePLM='plm-prim'},    	-- 0.00061365378346838
			--{solver='euler-hllc', hllcMethod=0, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons-alone'},	-- 0.00061365378346796
			--{solver='euler-hllc', hllcMethod=0, integrator='Runge-Kutta 3, TVD'},							-- 0.00061073285983012
			--{solver='euler-hllc', hllcMethod=0, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons'},		-- 0.00061073285983012

			-- hllcMethod == 1
			--{solver='euler-hllc', hllcMethod=1, integrator='forward Euler', usePLM='plm-cons-alone'},   -- 0.00060445812904152
			--{solver='euler-hllc', hllcMethod=1, integrator='forward Euler', usePLM='plm-prim'},         -- 0.00060393770417544
			--{solver='euler-hllc', hllcMethod=1, integrator='forward Euler'},							-- 0.00048873499978869
			--{solver='euler-hllc', hllcMethod=1, integrator='forward Euler', usePLM='plm-cons'},         -- 0.00048873499978869
			--{solver='euler-hllc', hllcMethod=1, integrator='Runge-Kutta 3, TVD', usePLM='plm-prim'},    -- 0.0006136537834679
			--{solver='euler-hllc', hllcMethod=1, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons-alone'}--, 0.00061365378346811
			--{solver='euler-hllc', hllcMethod=1, integrator='Runge-Kutta 3, TVD'},						-- 0.00061073285983055
			--{solver='euler-hllc', hllcMethod=1, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons'},	-- 0.00061073285983055

			-- hllcMethod == 2 (default)
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-cons-alone'},     	-- 0.00060447897937016
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-prim'},           	-- 0.00060393770371531
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler'},								-- 0.00048873499978778
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-cons'},           	-- 0.00048873499978778
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-prim'},      	-- 0.0006136537834683
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons-alone'},	-- 0.00061365378346809	
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD'},							-- 0.00061073285982937
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons'},	    -- 0.00061073285982937

			-- hllcMethod == 2 (default), but only using left and right state eigenvalues for waves -- not interface at all
			--  looks the same as the interface-based method, so ... why even use the interface waves?
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-cons-alone'},     	-- 0.00060449481040124
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-prim'},           	-- 0.00060393770360834
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler'},								-- 0.00048873499978885
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-cons'},           	-- 0.00048873499978885
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-prim'},      	-- 0.00061365378346831
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons-alone'},	-- 0.00061365378346856
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD'},							-- 0.0006107328598298
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons'},	    -- 0.0006107328598298

			-- hllcMethod == 2 (default), but without that last needless condition
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-cons-alone'},     	-- 0.00060487383014966 -- .1% worse, rest are the same
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-prim'},           	-- 0.00060393770418051
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler'},									-- 0.00048873499978941
			--{solver='euler-hllc', hllcMethod=2, integrator='forward Euler', usePLM='plm-cons'},           	-- 0.00048873499978941
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-prim'},      	-- 0.00061365378346782
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons-alone'},	-- 0.00061365378346786
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD'},								-- 0.0006107328598304
			--{solver='euler-hllc', hllcMethod=2, integrator='Runge-Kutta 3, TVD', usePLM='plm-cons'},	    	-- 0.0006107328598304


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
			--{solver='roe', integrator='forward Euler', fluxLimiter='Lax-Wendroff'},		-- 1.2048891136515e-06

			-- my PLM attempts:
			--{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},			-- 0.00049148119638364
			--{solver='roe', integrator='forward Euler', usePLM='plm-eig'},						-- 0.00049103769947135
			--{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim'},				-- 0.00049075156058252
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons'},					-- 0.00048873499978677
			--{solver='roe', integrator='forward Euler', usePLM='plm-athena'},					-- 0.00014403561237557

			-- various explicit integrators:
			--  once again, not so great
			--{solver='roe', integrator='Runge-Kutta 3', fluxLimiter='Lax-Wendroff'},			-- 0.063626499155696
			--{solver='roe', integrator='Runge-Kutta 3, TVD', fluxLimiter='Lax-Wendroff'},		-- 0.00012292610695107
			--{solver='roe', integrator='Runge-Kutta 4, non-TVD', fluxLimiter='Lax-Wendroff'},	-- 0.0001229260915465
			--{solver='roe', integrator='Runge-Kutta 4, TVD', fluxLimiter='Lax-Wendroff'},		-- 0.00012292609150301
			--{solver='roe', integrator='Runge-Kutta 4', fluxLimiter='Lax-Wendroff'},			-- 0.0001229260914525
			--{solver='roe', integrator='Runge-Kutta 2, non-TVD', fluxLimiter='Lax-Wendroff'},	-- 0.00012292557063657
			--{solver='roe', integrator='Runge-Kutta 2, TVD', fluxLimiter='Lax-Wendroff'},		-- 0.00012292557063653
			--{solver='roe', integrator='Runge-Kutta 2 Heun', fluxLimiter='Lax-Wendroff'},		-- 0.00012292557063543
			--{solver='roe', integrator='Runge-Kutta 2 Ralston', fluxLimiter='Lax-Wendroff'},	-- 0.00012292557063505
			--{solver='roe', integrator='Runge-Kutta 2', fluxLimiter='Lax-Wendroff'},			-- 0.00012292557063748
			--{solver='roe', integrator='Runge-Kutta 4, 3/8ths rule', fluxLimiter='Lax-Wendroff'},-- 3.0773114456088e-05
			
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
		local t = solver.t
		return xs:map(function(x,i)
			return rho0 + rho1 * math.sin(k0 * (x - u0 * t))
		end)
	end,

	-- TODO make sure solver_t->init_v0x == 1/duration and solver_t->maxs.x - mins.x == 2
	-- otherwise, for durations t=100 and t=1 the results look close enough to the same
	-- or just use what the Mara demo had:
	duration = 0.1,
	
	mins = {0,0,0},
	maxs = {1,1,1},
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
		{	-- (these numbers are all for duration=.1)
			-- schemes
			--{solver='weno5', integrator='forward Euler'},									-- 0.028050334485117
			--{solver='hll', integrator='forward Euler'},									-- 0.0039886633966807
			--{solver='euler-hllc', integrator='forward Euler'},							-- 0.0036984733332097
			--{solver='euler-burgers', integrator='forward Euler'},							-- 0.0031694650615551

			-- flux-limiters w/roe scheme:
			--{solver='roe', integrator='forward Euler', fluxLimiter='Lax-Wendroff'},		-- 0.47968895872093		-- the most accurate of the sine test, the least accurate of the Sod test ...
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
			
			-- plm-cons / slopeLimiter has funny behavior with Sod: it dips down but then gets worse.  
			-- I wouldn't be surprised if my exact values are off, so maybe that's why
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='CHARM'},				-- 0.55596019731028
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Lax-Wendroff'},		-- 0.47968895872093
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Beam-Warming'},		-- 0.47968895872093
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='van Leer'},			-- 0.47968895872093
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Fromm'},				-- 0.47968895872093
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Barth-Jespersen'},	-- 0.47968895872093
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='smart'},				-- 0.47968895872093
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='superbee'},			-- 0.010762308094614
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='monotized central'},	-- 0.007499966899957
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Koren'},				-- 0.006607365540715
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='HQUICK'},				-- 0.0064971868731382
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='HCUS'},				-- 0.0062065447431732
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='donor cell'},			-- 0.0056103088682637
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='UMIST'},				-- 0.0055130963289863
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Sweby'},				-- 0.0049485758799912
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='ospre'},				-- 0.0046940581345551
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='Oshker'},				-- 0.0039681296941796
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='van Albada 1'},		-- 0.0038128702530365
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='van Albada 2'},		-- 0.0037230670653778
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons', slopeLimiter='minmod'},				-- 0.0036575734195358

			-- other PLM attempts:
			--{solver='roe', integrator='forward Euler', usePLM='plm-cons'},				-- 0.0036660453357688
			--{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim-ref'},		-- 0.0023632029528882	-- plm-eig-prim-ref is worse than plm-eig-prim
			--{solver='roe', integrator='forward Euler', usePLM='plm-eig-prim'},			-- 0.0021682982996781	-- plm-eig-prim is worse than plm-eig
			--{solver='roe', integrator='forward Euler', usePLM='plm-eig'},					-- 0.0019239440708115
			--{solver='roe', integrator='forward Euler', usePLM='plm-athena'},				-- 0.0016947738688964
			--{solver='roe', integrator='forward Euler', usePLM='ppm-experimental'},		
			
			-- various explicit integrators:
			-- RK2 is exactly the same as forward Euler ... which means it's broken
			-- and everything else is worse														  before bug fix		after
			--{solver='roe', integrator='Runge-Kutta 2, non-TVD', fluxLimiter='superbee'},		-- 0.30989158592001		-- 0.47968895872093
			--{solver='roe', integrator='Runge-Kutta 3', fluxLimiter='superbee'},				-- 0.028049852882381    -- 0.038331373202177
			--{solver='roe', integrator='Runge-Kutta 2', fluxLimiter='superbee'},				-- 0.00187534979448	    -- 0.0036943914318805
			--{solver='roe', integrator='Runge-Kutta 2 Ralston', fluxLimiter='superbee'},		-- 0.014905909835352    -- 0.0036934820748332
			--{solver='roe', integrator='Runge-Kutta 2, TVD', fluxLimiter='superbee'},			-- 0.0277197399486      -- 0.0036928026259879
			--{solver='roe', integrator='Runge-Kutta 2 Heun', fluxLimiter='superbee'},			-- 0.0277197399486      -- 0.0036928026259879
			--{solver='roe', integrator='Runge-Kutta 4', fluxLimiter='superbee'},				-- 0.010605564924952    -- 0.0036907642308442
			--{solver='roe', integrator='Runge-Kutta 4, non-TVD', fluxLimiter='superbee'},		-- 0.010605564924952    -- 0.0036907642308437
			--{solver='roe', integrator='Runge-Kutta 4, TVD', fluxLimiter='superbee'},			-- 0.010605515773635    -- 0.0036907403381736
			--{solver='roe', integrator='Runge-Kutta 3, TVD', fluxLimiter='superbee'},			-- 0.019060751585123    -- 0.0036901766164101
			--{solver='roe', integrator='Runge-Kutta 4, 3/8ths rule', fluxLimiter='superbee'},	-- 0.0083893332244032   -- 0.0034902119354803
	
			-- various RK's with PLM
			-- even though plm-athena is the best among forward-Euler,
			--  it doesn't seem to do too well with various RK integrators
			--{solver='roe', integrator='Runge-Kutta 2, non-TVD', usePLM='plm-athena'},		-- 0.50027377262497
			--{solver='roe', integrator='Runge-Kutta 3', usePLM='plm-athena'},				-- 0.038327167106819
			--{solver='roe', integrator='Runge-Kutta 2', usePLM='plm-athena'},				-- 0.0035544851041015
			--{solver='roe', integrator='Runge-Kutta 2 Ralston', usePLM='plm-athena'},		-- 0.0035579071138359
			--{solver='roe', integrator='Runge-Kutta 2, TVD', usePLM='plm-athena'},			-- 0.0035580548511293
			--{solver='roe', integrator='Runge-Kutta 2 Heun', usePLM='plm-athena'},			-- 0.0035580548511293
			--{solver='roe', integrator='Runge-Kutta 4', usePLM='plm-athena'},				-- 0.0035465402301359
			--{solver='roe', integrator='Runge-Kutta 4, non-TVD', usePLM='plm-athena'},		-- 0.0035465402301358
			--{solver='roe', integrator='Runge-Kutta 4, TVD', usePLM='plm-athena'},			-- 0.0035465949424922
			--{solver='roe', integrator='Runge-Kutta 3, TVD', usePLM='plm-athena'},			-- 0.0035456783040898
			--{solver='roe', integrator='Runge-Kutta 4, 3/8ths rule', usePLM='plm-athena'},	-- 0.0035367821065472

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

	mins = {-1,-1,-1},
	maxs = {1,1,1},
	duration = .2,
}

local function calcError(exact, ys)
	local diff = matrix(exact) - matrix(ys)
	local err = diff:normL1() / #diff
	return err
end

local problem = problems[problemName]

local testdatas = table()
local errorsForConfig = table()
local errorNames = table()

-- for history and compare, this is the size used
--local singleSize = 1024
local singleSize = 64

local dim = 1
local sizes = (plotCompare or plotErrorHistory) and table{singleSize} or range(3,10):map(function(x) return 2^x end)
for _,cfg in ipairs(problem.configurations) do
	cfg = table(cfg)

	local destName = string.trim(tolua(cfg):match('^{(.*)}$'):gsub('%s+', ' '):gsub('"', ''))
print(destName)

	local destFilename = destName
		:gsub('/', '')
		:gsub('{', '(')
		:gsub('}', ')')
	
	cfg.dim = dim
	cfg.cfl = .6
	cfg.coord = 'cartesian'
	cfg.mins = problem.mins
	cfg.maxs = problem.maxs

	--[[
	data cached per-test:
	size[size]
		.xs[index] 
		.ys[index] 
	--]]
	local testdata
	local srcfn = rundir..'/'..destFilename..'.lua'
	local srcfiledata = file[srcfn]
	if srcfiledata then
		testdata = fromlua(srcfiledata)
	end
	testdata = testdata or {}
	do	--if nocache or not testdata.size then
		testdata.size = testdata.size or table()
		for _,size in ipairs(sizes) do
			
			testdata.size[size] = testdata.size[size] or table()
					
			if nocache 
			or not testdata.size[size].xs
			or not testdata.size[size].ys
			-- or either no exact or uselin
			then
				cfg = table(cfg)
				cfg.gridSize = {size}
print()
print(size)
print()	

				if plotErrorHistory then
					testdata.size[size].ts = table()
					testdata.size[size].errorsForTime = table()
				end

				local duration = tonumber(problem.duration) or error("expected problem.duration")
				
				local function getSolverGraph(solver)
					local _, var = solver.displayVars:find(nil, function(var) return var.name == 'U rho' end)
					assert(var, "failed to find U rho var")
					solver:calcDisplayVarToBuffer(var)	
					-- now in solver.reduceBuf
					local ptr = solver.calcDisplayVarToTexPtr
					local numCells = solver.numCells
					local numGhost = solver.numGhost
					local app = solver.app
					app.cmds:enqueueReadBuffer{buffer=solver.reduceBuf, block=true, size=ffi.sizeof(app.real) * numCells, ptr=ptr}
					-- now in ptr
					
					local xs = range(numCells-2*numGhost):mapi(function(i)
						return (i-.5) * solver.solverPtr.grid_dx.x + solver.solverPtr.mins.x
					end)
					local ys = range(numCells-2*numGhost):mapi(function(i)
						return ptr[i+numGhost-1]
					end)
					local exact
					if not uselin then
						exact = problem:makeExact(xs, solver)
					end			
					return xs, ys, exact
				end

				local function calcAndSaveTimeAndError(solver)
					local xs, ys, exact = getSolverGraph(solver)
					testdata.size[size].ts:insert(solver.t)
					local err = calcError(exact, ys)
					testdata.size[size].errorsForTime:insert(err)
				end
				
				local App = class(require 'app')
				function App:setup(clArgs)
					cfg.app = self
					local solver = require('solver.'..cfg.solver)(cfg)
					self.solvers:insert(solver)
					self.exitTime = duration
					self.running = true
					if plotErrorHistory then
						local oldupdate = solver.update
						solver.update = function(...)
							calcAndSaveTimeAndError(solver)
							return oldupdate(...)
						end
					end
				end
				
				function App:requestExit()
					App.super.requestExit(self)
				
					-- now compare the U buffer to the exact 
					assert(#self.solvers == 1)
					local solver = self.solvers[1]
					local xs,ys,exact = getSolverGraph(solver)
					testdata.size[size].xs = xs
					testdata.size[size].ys = ys
					testdata.size[size].exact = exact
					
					local err = calcError(exact, ys)
					testdata.size[size].error = err
					if plotErrorHistory then
						testdata.size[size].ts:insert(solver.t)
						testdata.size[size].errorsForTime:insert(err)
					end
				end	
				
				local app  = App()
				app:run()
			end
			
			local xs = setmetatable(assert(testdata.size[size].xs), table)
			local ys = setmetatable(assert(testdata.size[size].ys), table)
			local exact = not uselin and assert(testdata.size[size].exact)
			
			if uselin then
				-- just use log/log regression to estimate where the best would be
				-- technically this won't be best, because most our samples tend to flatten at the bottom
				-- and in that case, the regression will point to a less accurate place than where it would if a smaller subset was used 
				local xavg = xs:sum() / #xs
				local yavg = ys:sum() / #ys

				local b1 = range(#xs):map(function(i)
					return (xs[i] - xavg) * (ys[i] - yavg)
				end):sum() / range(#xs):map(function(i)
					return (xs[i] - xavg)^2
				end):sum()
				local b0 = yavg - b1 * xavg
				exact = xs:map(function(x)
					return b0 + b1 * x
				end)
			end

			if plotCompare then -- plotting immediately
				gnuplot{
					output = rundir..'/compare-graphs.png',
					style = 'data lines',
					data = {xs, ys, exact},
					{using = '1:2', title=''..size},
					exact and {using = '1:3', title='exact'} or nil,
				}
				os.exit()
			end

			if file.stop then file.stop = nil os.exit(1) end
		end
	end
	testdata.name = destName
	file[srcfn] = tolua(testdata)
local errors = sizes:map(function(size) return testdata.size[size].error end)
print(table.last(errors))
	errorsForConfig:insert(errors)
	errorNames:insert(destName)
	testdatas:insert(testdata) 
end

-- [[ plot errors per size
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

-- [[ plot error histories
if plotErrorHistory then
	local size = singleSize 
	local data = table()
	for _,testdata in ipairs(testdatas) do
		data:insert(assert(testdata.size[size].ts))
		data:insert(assert(testdata.size[size].errorsForTime))
	end
	gnuplot(table(
		{
			output = rundir..'/error-history.png',
			style = 'data lines',
			log = 'xy',
			data = data,
		},
		testdatas:map(function(testdata,i)
			return {using=(2*i-1)..':'..(2*i), title=testdata.name}
		end)
	))
end
--]]


local dim = cmdline.dim or 1
local args = {
	app = self, 
	eqn = cmdline.eqn,
	dim = dim,
	
	integrator = cmdline.integrator or 'forward Euler',	
	--integrator = 'Iterative Crank-Nicolson',
	--integrator = 'Runge-Kutta 2',
	--integrator = 'Runge-Kutta 2 Heun',
	--integrator = 'Runge-Kutta 2 Ralston',
	--integrator = 'Runge-Kutta 3',
	--integrator = 'Runge-Kutta 4',
	--integrator = 'Runge-Kutta 4, 3/8ths rule',
	--integrator = 'Runge-Kutta 2, TVD',
	--integrator = 'Runge-Kutta 2, non-TVD',
	--integrator = 'Runge-Kutta 3, TVD',
	--integrator = 'Runge-Kutta 4, TVD',
	--integrator = 'Runge-Kutta 4, non-TVD',
	--integrator = 'backward Euler',	-- The epsilon on this is very sensitive.  Too small and it never converges.  Too large and it stops convergence too soon.
	--integratorArgs = {verbose=true},

	--fixedDT = .0001,
	cfl = cmdline.cfl or .6/dim,	-- 1/dim,
	
	fluxLimiter = cmdline.fluxLimiter or 'superbee',
	--fluxLimiter = 'monotized central',
	--fluxLimiter = 'donor cell',
	
	-- piecewise-linear slope limiter
	-- TODO rename this to 'calcLR' or something
	--usePLM = 'plm-cons',
	--usePLM = 'plm-cons-alone',
	--usePLM = 'plm-prim-alone',
	--usePLM = 'plm-eig',
	--usePLM = 'plm-eig-prim',
	--usePLM = 'plm-eig-prim-ref',
	--usePLM = 'plm-athena',			-- based on Athena.  most accurate from 1D sod tests atm
	--usePLM = 'ppm-experimental',	-- FIXME one more attempt to figure out all the PLM stuff, based on 2017 Zingale
	--usePLM = 'weno',				-- TODO make WENO one of these 'usePLM' methods. rename it to 'construct LR state method' or something.  then we can use CTU with WENO.
	
	-- only enabled for certain usePLM methods
	--slopeLimiter = 'minmod',

	-- this is functional without usePLM, but doing so falls back on the cell-centered buffer, which with the current useCTU code will update the same cell twice from different threads
	--useCTU = true,
	
	--[[ Cartesian
	coord = 'cartesian',
	mins = cmdline.mins or {-1, -1, -1},
	maxs = cmdline.maxs or {1, 1, 1},

	-- 256^2 = 2^16 = 2 * 32^3
	gridSize = (
		({ 	-- size options based on OpenCL vendor ...
			['NVIDIA CUDA/GeForce GTX 1080 Ti'] = {
				{256,1,1},
				{256,256,1},
				{64,64,64},
			},
			['NVIDIA CUDA/GeForce GTX 1080'] = {
				{256,1,1},
				{256,256,1},
				{32,32,32},
			},
			['Intel(R) OpenCL/Intel(R) HD Graphics 520'] = {
				{256,1,1},
				{64,64,1},
				{16,16,16},
			},
			['Intel(R) OpenCL HD Graphics/Intel(R) Gen9 HD Graphics NEO'] = {
				{64,1,1},
				{64,64,1},
				{16,16,16},
			},
		})[platformName..'/'..deviceName] 
		-- default size options
		or {
			{256,1,1},
			{256,256,1},
			{32,32,32},
		}
	)[dim],
	boundary = {
		xmin=cmdline.boundary or 'freeflow',
		xmax=cmdline.boundary or 'freeflow',
		ymin=cmdline.boundary or 'freeflow',
		ymax=cmdline.boundary or 'freeflow',
		zmin=cmdline.boundary or 'freeflow',
		zmax=cmdline.boundary or 'freeflow',
	},
	--]]
	--[[ cylinder
	coord = 'cylinder',
	coordArgs = {anholonomic=true},			-- disable to use non-physical, holonomic coordinates
	mins = cmdline.mins or {.2, 0, -.25},
	maxs = cmdline.maxs or {1, 2*math.pi, .25},
	gridSize = ({
		{128, 1, 1}, -- 1D
		{64, 256, 1}, -- 2D
		{32, 32, 32}, -- 3D
	})[dim],
	boundary = {
		-- r
		xmin=cmdline.boundary or 'mirror',		-- hmm, how to treat the r=0 boundary ...
		xmax=cmdline.boundary or 'freeflow',
		-- theta
		ymin=cmdline.boundary or 'periodic',
		ymax=cmdline.boundary or 'periodic',
		-- z
		zmin=cmdline.boundary or 'freeflow',
		zmax=cmdline.boundary or 'freeflow',
	},
	--]]
	-- [[ Sphere: r, θ, φ 
	coord = 'sphere',
	--coordArgs = {volumeDim = 3},	-- use higher dimension volume, even if the grid is only 1D to 3D
	mins = cmdline.mins or {1, 0, -math.pi},
	maxs = cmdline.maxs or {10, math.pi, math.pi},
	gridSize = ({
		{256, 1, 1}, -- 1D
		{32, 32, 1}, -- 2D
		{16, 16, 16}, -- 3D
	})[dim],
	boundary = {
		xmin=cmdline.boundary or 'freeflow',	-- how to handle?
		xmax=cmdline.boundary or 'freeflow',
		ymin=cmdline.boundary or 'freeflow',	-- ?
		ymax=cmdline.boundary or 'freeflow',	-- ?
		zmin=cmdline.boundary or 'periodic',
		zmax=cmdline.boundary or 'periodic',
	},
	--]]
	-- symbolic math gets too complex
	-- this is a good reason for arbitrary mesh support
	--[[ torus
	coord = 'torus',
	-- hmm, right now sphere's variables change per-dimension used ...
	mins = cmdline.mins or {0, 0, 0},
	maxs = cmdline.maxs or {2*math.pi, 2*math.pi, 1},
	gridSize = {
		cmdline.gridSize or 16,
		cmdline.gridSize or 16,
		cmdline.gridSize or 16,
	},
	boundary = {
		xmin=cmdline.boundary or 'periodic',
		xmax=cmdline.boundary or 'periodic',
		ymin=cmdline.boundary or 'periodic',
		ymax=cmdline.boundary or 'periodic',
		zmin=cmdline.boundary or 'mirror',
		zmax=cmdline.boundary or 'mirror',
	},
	--]]

	--useGravity = true,

	-- TODO separate initStates for each class of equation
	-- this would cohese better with the combined solvers
	-- i.e. a fluid initState, a Maxwell init-state, and a GR init-state
	-- ... but that means splitting the MHD init-states across M and HD ...
	-- how about just stacking initStates?
	-- and letting each one assign what values it wants.
	-- still to solve -- how do we specify initStates for combined solvers with multiple sets of the same variables (ion/electron, etc)

	-- no initial state means use the first
	--initState = cmdline.initState,
	
	-- Euler / SRHD / MHD initial states:
	--initState = 'constant',
	--initState = 'constant with velocity',
	--initState = 'linear',
	--initState = 'gaussian',
	--initState = 'advect wave',
	--initState = 'sphere',
	--initState = 'rarefaction wave',
	
	--initState = 'Sod',
	--initState = 'rectangle',
	--initState = 'Sedov',
	--initState = 'Noh',
	--initState = 'implosion',
	--initState = 'Kelvin-Helmholtz',
	--initState = 'Rayleigh-Taylor',	--FIXME
	--initState = 'Colella-Woodward',
	--initState = 'double mach reflection',
	--initState = 'square cavity',
	--initState = 'shock bubble interaction',		-- with usePLM only works with prim or with athena
	--initState = 'Richmyer-Meshkov',

	--initState = 'configuration 1',
	--initState = 'configuration 2',
	--initState = 'configuration 3',
	--initState = 'configuration 4',
	--initState = 'configuration 5',
	--initState = 'configuration 6',

	-- self-gravitation tests:
	--initState = 'self-gravitation - Earth',	-- validating units along with self-gravitation.  TODO this doesn't produce 9.8 m/s^2 just yet
	--initState = 'self-gravitation test 1',
	--initState = 'self-gravitation test 1 spinning',
	--initState = 'self-gravitation test 2',		--FIXME
	--initState = 'self-gravitation test 2 orbiting',
	--initState = 'self-gravitation test 4',
	--initState = 'self-gravitation soup',	--FIXME
	
	-- those designed for SRHD / GRHD:
	--initState = 'relativistic shock reflection',			-- FIXME.  these initial conditions are constant =P
	--initState = 'relativistic blast wave test problem 1',
	--initState = 'relativistic blast wave test problem 2',
	--initState = 'relativistic blast wave interaction',		-- in 2D this only works with no limiter / lots of dissipation 

	-- states for ideal MHD or two-fluid (not two-fluid-separate)
	--initState = 'Brio-Wu',
	--initState = 'Orszag-Tang',
	--initState = 'MHD rotor',
	--initState = 'GEM challenge', eqnArgs = {useEulerInitState=false},
	--initState = 'spinning magnetic fluid',
	--initState = 'magnetic fluid',
	--initState = '2017 Degris et al',
	--initState = 'that one mhd simulation from youtube',
	
	-- 2002 Dedner
	--initState = '2002 Dedner peak Bx',
	--initState = '2002 Dedner 1D Riemann',
	--initState = '2002 Dedner Shock Reflection',
	--initState = '2002 Dedner 2D Riemann problem',
	--initState = '2002 Dedner Kelvin-Helmholtz',

	-- Maxwell:
	--initState = 'Maxwell default',
	--initState = 'Maxwell empty waves',
	--initState = 'Maxwell scattering around cylinder',
	--initState = 'Maxwell scattering around pyramid',
	--initState = 'Maxwell scattering around square',
	--initState = 'Maxwell scattering around Koch snowflake',
	--initState = 'Maxwell wire',
	--initState = 'Maxwell transverse waves',
	
	-- hmm, I think I need a fluid solver for this, not just a Maxwell solver ...
	--initState = 'Maxwell Lichtenberg',	

	-- Maxwell+HD
	--initState = 'two-fluid emhd modified Brio-Wu',
	--initState = 'two-fluid EMHD soliton ion',
	--initState = 'two-fluid EMHD soliton electron',
	--initState = 'two-fluid EMHD soliton maxwell',


	-- Einstein
	initState = 'Minkowski',
	--initState = 'gaussian perturbation',
	--initState = 'plane gauge wave',


	--initState = 'Alcubierre warp bubble',
	
	--initStateArgs = {R=.5, sigma=8, speed=.1},	-- sub-luminal
	
	--initStateArgs = {R=.5, sigma=8, speed=1.1},		-- super-luminal 1.1x
	-- ... works with
	--	size=64x64 solver=adm3d int=fe plm=athena ctu
	--  size=64x64 solver=adm3d int=fe plm=athena
	--  size=64x64 solver=adm3d int=fe flux-limiter=superbee
	
	--initStateArgs = {R=.5, sigma=8, speed=2},		-- super-luminal 2x
	-- ... works with
	--  size=64x64 solver=adm3d int=fe flux-limiter=superbee
	
	--initStateArgs = {R=.5, sigma=8, speed=10},		-- super-luminal 10x
	--  size=64x64 solver=roe eqn=adm3d int=fe flux-limiter=superbee ... eventually explodes
	--  size=64x64 solver=roe eqn=bssnok int=be flux-limiter=superbee ... eventually explodes as well
	--  size=128x128 solver=hll eqn=adm3d int=fe flux-limiter=superbee ... runs for a really long time 

	
	--initState = 'black hole - Schwarzschild',
	
	
	--initState = 'black hole - isotropic',	-- this one has momentum and rotation and almost done with multiple sources.  TODO parameterize
	
	--initState = 'black hole - SENR/NumPy',

	--[[ single black hole, stationary
	initStateArgs = {
		bodies = {
			{
				R = .1,
				P_u = {0,0,0},
				S_u = {0,0,0},
				pos = {0,0,0},
			},
		},
	},
	--]]

	--[[ single black hole, spinning, demonstrating ergosphere formation
	initStateArgs = {
		bodies = {
			{
				R = .0001,
				P_u = {0,0,0},
				S_u = {0,0,1},
				pos = {0,0,0},
			},
		},
	},
	--]]
	--[[ single black hole 
	-- from 2006 Brugmann et al - "Calibration of a Moving Puncture Simulation" (though I'm not using moving puncture method)
	-- set to 1/10th the scale 
	initStateArgs = {
		bodies = {
			{
				R = .0505,
				P_u = {0,0,0},
				S_u = {0,0,.0866},
				pos = {0,0,0},
			},
		},
	},
	--]]

	--[[ single black hole, boosted
	initStateArgs = {
		bodies = {
			{
				R = .05,
				P_u = {.01,0,0},
				S_u = {0,0,0},
				pos = {0,0,0},
			},
		},
	},
	--]]

	--[[ single black hole, spinning and boosted perpendicular to axis
	initStateArgs = {
		bodies = {
			{
				R = .05,
				P_u = {.01,0,0},
				S_u = {0,0,1},
				pos = {0,0,0},
			},
		},
	},
	--]]

	--[[ single black hole, spinning and boosted along axis
	initStateArgs = {
		bodies = {
			{
				R = .05,
				P_u = {0,0,.01},
				S_u = {0,0,1},
				pos = {0,0,0},
			},
		},
	},
	--]]

	--[[ binary black hole, head-on collision 
	initStateArgs = {
		bodies = {
			{
				R = .01,
				P_u = {-.01,0,0},
				S_u = {0,0,0},
				pos = {.5,0,0},
			},
			{
				R = .01,
				P_u = {.01,0,0},
				S_u = {0,0,0},
				pos = {-.5,0,0},
			},
		},
	},
	--]]

	--[[ binary black hole from 2006 Brugmann et al, same as above
	initStateArgs = {
		bodies = {
			{
				R = .0505,
				P_u = {0,.0133,0},
				S_u = {0,0,.0866},
				pos = {.3257,0,0},
			},
			{
				R = .0505,
				P_u = {0,-.0133,0},
				S_u = {0,0,.0866},
				pos = {-.3257,0,0},
			},
		},
	},
	--]]


	--initState = 'stellar model',
	--initState = '1D black hole - wormhole form',

	
	--initState = 'Gowdy waves',
	--initState = 'testbed - robust',	-- not working with fv solvers
	--initState = 'testbed - gauge wave',	-- not working with forward-euler finite-difference solvers
	--initState = 'testbed - gauge wave - diagonal',
	--initState = 'testbed - linear wave',
	--initState = 'testbed - linear wave - diagonal',
	--initState = 'testbed - Gowdy',


	-- NLS
	--initState = 'Gaussian',
	--initState = 'Ring',
	--initState = 'Oscillatory',
}


if cmdline.solver then self.solvers:insert(require('solver.'..cmdline.solver)(table(args, cmdline))) return end


-- wave equation


--self.solvers:insert(require 'solver.roe'(table(args, {eqn='wave'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='wave'})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='wave', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='wave', wenoMethod='2008 Borges', order=5})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='wave', wenoMethod='2010 Shen Zha', order=5})))


-- compressible Euler equations


--self.solvers:insert(require 'solver.roe'(table(args, {eqn='euler'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='euler'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='euler'})))

--self.solvers:insert(require 'solver.euler-hllc'(table(args, {hllcMethod=0})))
--self.solvers:insert(require 'solver.euler-hllc'(table(args, {hllcMethod=1})))
--self.solvers:insert(require 'solver.euler-hllc'(table(args, {hllcMethod=2})))

-- NOTICE, these are very accurate with RK4, etc., but incur oscillations with Forward-Euler
-- TODO weno doesn't seem to work with self-gravitation
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='2008 Borges', order=5})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=5})))

--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=7})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='2008 Borges', order=7})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=7})))

-- order=9 and above 2D Sod test, things tend to explode ... maybe I should take away the betaCoeffs denominator (since they get normalized anyways?)
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=9})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='2008 Borges', order=9})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=9})))

--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=11})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='2008 Borges', order=11})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=11})))

-- hmm, 2D Sod 64x64 RK4 fails at just past 1 second ... 
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=13})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='2008 Borges', order=13})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=13})))

-- blows up.  maybe I need an implicit RK scheme...
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=13, integrator='backward Euler'})))

-- still haven't added source terms to this
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='navstokes-wilcox'})))


-- compressible Euler equations - Burgers solver


-- f.e. and b.e. are working, but none of the r.k. integrators
-- PLM isn't implemented yet
-- neither is source term / poisson stuff
--self.solvers:insert(require 'solver.euler-burgers'(args))


-- special relativistic compressible hydrodynamics


-- rel blast wave 1 & 2 works in 1D at 256 with superbee flux lim
-- rel blast wave interaction works with superbee flux lim in 1D works at 256, fails at 1024 with float (works with double)
-- 	256x256 double fails with F.E., RK2-Heun, RK2-Ralston, RK2-TVD, RK3, RK4-3/8ths,
-- rel blast wave 1 doesn't work in 64x64. with superbee flux lim
-- rel blast wave 2 with superbee flux lim, Roe solver, works at 64x64 with forward euler
-- 	at 256x256 fails with F.E, RK2, RK2-non-TVD., RK3-TVD, RK4, RK4-TVD, RK4-non-TVD 
--    but works with RK2-Heun, RK2-Ralston, RK2-TVD, RK3, RK4-3/8ths
-- Kelvin-Helmholtz works for all borderes freeflow, float precision, 256x256, superbee flux limiter
--self.solvers:insert(require 'solver.srhd-roe'(args))
	-- TODO can't use these until I get eigen_forInterface working in eqn/srhd.cl
--self.solvers:insert(require 'solver.srhd-hll'(args))
--self.solvers:insert(require 'solver.srhd-weno'(table(args, {eqn='srhd', wenoMethod='2010 Shen Zha', order=5})))


-- general relativistic compressible hydrodynamics


-- TODO remove calcEigenBasis from eqn/grhd.cl
-- this is the solver with plug-ins for ADM metric, 
-- yet doesn't come coupled with any other solver, so it will act just like a srhd solver
--self.solvers:insert(require 'solver.grhd-roe'(args))


-- ideal magnetohydrodynamics


-- with superbee flux lim:  
-- Brio-Wu works in 1D at 256, works in 2D at 64x64 in a 1D profile in the x and y directions.
-- Orszag-Tang with forward Euler integrator fails at 64x64 around .7 or .8
-- 		but works with 'Runge-Kutta 4, TVD' integrator at 64x64
-- 		RK4-TVD fails at 256x256 at just after t=.5
--		and works fine with backwards Euler 
-- when run alongside HD Roe solver, curves don't match (different heat capacity ratios?)
--		but that could be because of issues with simultaneous solvers.
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='mhd'})))

-- this runs, but of course it's missing a few waves ...
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='mhd'})))

-- explodes with Orszag-Tang
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='mhd', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='mhd', wenoMethod='2010 Shen Zha', order=5})))

--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='mhd'})))
-- TODO FIXME
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='glm-mhd'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='glm-mhd'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='glm-mhd'})))


-- Maxwell


--self.solvers:insert(require 'solver.roe'(table(args, {eqn='maxwell'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='maxwell'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='maxwell'})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='maxwell', wenoMethod='1996 Jiang Shu', order=5})))


-- GLM Maxwell


--self.solvers:insert(require 'solver.roe'(table(args, {eqn='glm-maxwell'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='glm-maxwell'})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='glm-maxwell', wenoMethod='2010 Shen Zha', order=7})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='glm-maxwell', wenoMethod='2010 Shen Zha', order=13})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='glm-maxwell'})))


-- Maxwell+HD two-fluid electron/ion solver
-- TODO FIXME
-- TODO, with the separate solver, use hll, so the ion, electron, and maxwell all use hll separately
-- TODO I made it even more difficult to implement with the addition of these real_ and cplx_ macros...
--self.solvers:insert(require 'solver.twofluid-emhd-separate-roe'(args))

-- ...so to try and get around that, here the two are combined into one solver:
-- it is stable, however since all variables are tied together, it integrates them together, which means everything is explicit-integration updated
-- ...which means, with the Maxwell equations waves propagating at the speed of light, that it goes very slow
-- TODO: I suppose I could make this work with my integrator by (1) removing the maxwell terms from the integration variable list and (2) providing a separate operator that updates them implicitly
-- TODO still needs PLM support
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='twofluid-emhd'})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='twofluid-emhd', wenoMethod='1996 Jiang Shu', order=9})))	-- exploded...
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='twofluid-emhd', wenoMethod='2010 Shen Zha', order=5})))


-- here's another one: two-fluid emhd with de Donder gauge linearized general relativity
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='twofluid-emhd-lingr'})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='twofluid-emhd-lingr', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='twofluid-emhd-lingr', wenoMethod='2010 Shen Zha', order=7})))


-- general relativity


--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm1d_v1'})))
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm1d_v2'})))
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d'})))
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d', eqnArgs={noZeroRowsInFlux=false}})))
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d', eqnArgs={useShift='MinimalDistortionElliptic'}})))	-- TODO finish me
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d', eqnArgs={useShift='MinimalDistortionEllipticEvolve'}})))	-- TODO finish me
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d', eqnArgs={useShift='2005 Bona / 2008 Yano'}})))	-- TODO finish me
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d', eqnArgs={useShift='HarmonicShiftCondition-FiniteDifference'}})))	-- breaks, even with b.e. integrator
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d', eqnArgs={useShift='LagrangianCoordinates'}})))	-- TODO finish me
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='z4_2008yano'})))
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='z4'})))

--self.solvers:insert(require 'solver.hll'(table(args, {eqn='adm1d_v1'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='adm1d_v2'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='adm3d'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='adm3d', eqnArgs={noZeroRowsInFlux=false}})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='adm3d', eqnArgs={useShift='HarmonicShiftCondition-FiniteDifference'}})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='z4'})))

--self.solvers:insert(require 'solver.weno'(table(args, {eqn='adm3d', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'solver.weno'(table(args, {eqn='adm3d', wenoMethod='2010 Shen Zha', order=7})))

--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='adm1d_v1', integrator='backward Euler'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='adm1d_v2', integrator='backward Euler'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='adm3d', integrator='backward Euler'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='z4_2008yano', integrator='backward Euler'})))



-- bssnok is working in 1D-3D Cartesian for RK4 but diverging for non-Cartesian
self.solvers:insert(require 'solver.bssnok-fd'(args))	-- default shift is HyperbolicGammaDriver
--self.solvers:insert(require 'solver.bssnok-fd'(table(args, {eqnArgs={useShift='none'}})))
--self.solvers:insert(require 'solver.bssnok-fd'(table(args, {eqnArgs={useShift='GammaDriver'}})))

-- Z4c finite difference, combining BSSNOK and Z4
-- FIXME something is asymmetric.  watch Theta.  Run warp bubble.
--self.solvers:insert(require 'solver.z4c-fd'(args))



-- GRHD+GR
-- here's the GRHD solver with the BSSNOK plugged into it
-- TODO inital condition root finding to make sure the EFE is satisifed
--self.solvers:insert(require 'solver.gr-hd-separate'(args))


-- GR+Maxwell.  params go to the Maxwell solver.

-- This is a Maxwell solver with the extra terms to account for GR.
-- It only simulates Maxwell -- it needs an additional GR solver to fully work correctly.
-- (Without a GR solver it will just operate in Minkowski spacetime.)
-- TODO this should be the einstein-maxwell eqns solver with flat space plugged into the metric
-- TODO rename to einstein-maxwell-roe
-- TODO it's unfinished
--self.solvers:insert(require 'solver.gr-maxwell-roe'(args))

-- TODO and this should be the einstein-maxwell with an einstein solver plugged in
-- TODO rename to einstein-maxwell-separate
--self.solvers:insert(require 'solver.gr-em-separate'(args))


-- nonlinear Schrodinger equation
--self.solvers:insert(require 'solver.nls'(args))


-- the start of unstructured meshes
--self.solvers:insert(require 'solver.meshsolver'(table(args, {eqn='euler', meshfile='n0012_113-33'})))
-- temp here -- to make sure ordinary solvers still run
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='euler', initState='Sod'})))


-- the start of AMR
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='euler'})))
--self.solvers:insert(require 'solver.amr'(require 'solver.roe')(table(args, {eqn='euler'})))

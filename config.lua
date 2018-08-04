local dim = 2
local args = {
	app = self, 
	eqn = cmdline.eqn,
	dim = cmdline.dim or dim,
	
	integrator = cmdline.integrator or 'forward Euler',	
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
	--integrator = 'backward Euler',
	
	--fixedDT = .0001,
	--cfl = .25/dim,
	
	fluxLimiter = cmdline.fluxLimiter or 'superbee',
	--fluxLimiter = 'monotized central',
	--fluxLimiter = 'donor cell',
	
	-- piecewise-linear slope limiter
	--usePLM = 'plm-cons',			-- works in conservative variable space, uses a slope limiter
	--usePLM = 'plm-eig',			-- works in conservative eigenspace, uses 2 slopes for the limiter (TODO incorporate slopeLimiter)
	--usePLM = 'plm-eig-prim',		-- works in primitive eigenspace, etc
	--usePLM = 'plm-eig-prim-ref',	-- works in primitive eigenspace, etc, subtracts out min & max.  doesn't work well with ideal mhd.
	--usePLM = 'plm-athena',		-- based on Athena, idk about this one
	--usePLM = 'ppm-experimental',	-- one more attempt to figure out all the PLM stuff, based on 2017 Zingale

	-- only enabled for certain usePLM methods
	--slopeLimiter = 'minmod',

	--useCTU = true,
	
	-- [[ Cartesian
	coord = 'cartesian',
	mins = cmdline.mins or {-1, -1, -1},
	maxs = cmdline.maxs or {1, 1, 1},

--[=[ for Richmyer-Meshkov
-- TODO let the init conds configure domain
mins = {-2,0,0},
maxs = {6,1,1},
--]=]

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
			['Intel(R) OpenCL/Intel(R) HD Graphics'] = {
				{256,1,1},
				{256,256,1},
				{32,32,32},
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
	-- TODO these next two seem very similar
	--[[ 1D radial
	coord = '1d_radial',
	mins = cmdline.mins or {0, 0, 0},
	maxs = cmdline.maxs or {1, 1, 1},
	gridSize = ({
		{256,1,1},
		{32,32,1},
		{32,32,32},
	})[dim],
	boundary = {
		xmin=cmdline.boundary or 'freeflow',
		xmax=cmdline.boundary or 'freeflow',
		ymin=cmdline.boundary or 'periodic',
		ymax=cmdline.boundary or 'periodic',
		zmin=cmdline.boundary or 'periodic',
		zmax=cmdline.boundary or 'periodic',
	},
	--]]
	--[[ sphere1d -- used for 1D radial profiles of spheres
	coord = 'sphere1d',
	mins = cmdline.mins or {.1, 0, -math.pi},
	maxs = cmdline.maxs or {1, math.pi, math.pi},
	gridSize = {
		cmdline.gridSize or 256,
		cmdline.gridSize or 128,
		cmdline.gridSize or 64,
	},
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
	-- for some reason [rmin, rmax] = [.5, 1] gets an explosion right at r=rmin, theta=0 ... but any other values work fine
	coord = 'cylinder',
	mins = cmdline.mins or {.1, 0, -.25},
	maxs = cmdline.maxs or {1, 2*math.pi, .25},
	gridSize = ({
		{128, 1, 1}, -- 1D
		{64, 256, 1}, -- 2D
		{16, 64, 16}, -- 3D
	})[dim],
	boundary = {
		-- r
		xmin=cmdline.boundary or 'mirror',		-- hmm, how to treat the r=0 boundary ...
		xmax=cmdline.boundary or 'mirror',
		-- theta
		ymin=cmdline.boundary or 'periodic',
		ymax=cmdline.boundary or 'periodic',
		-- z
		zmin=cmdline.boundary or 'mirror',
		zmax=cmdline.boundary or 'mirror',
	},
	--]]
	--[[ sphere
	coord = 'sphere',
	mins = cmdline.mins or {.1 * math.pi, -math.pi, .1},
	maxs = cmdline.maxs or {.9 * math.pi, math.pi, 1},
	gridSize = {
		cmdline.gridSize or 16,
		cmdline.gridSize or 16,
		cmdline.gridSize or 16,
	},
	boundary = {
		xmin=cmdline.boundary or 'mirror',
		xmax=cmdline.boundary or 'mirror',
		ymin=cmdline.boundary or 'periodic',
		ymax=cmdline.boundary or 'periodic',
		zmin=cmdline.boundary or 'mirror',
		zmax=cmdline.boundary or 'mirror',
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
	--initState = 'Sedov',
	--initState = 'Noh',
	--initState = 'implosion',
	--initState = 'Kelvin-Helmholtz',
	--initState = 'Rayleigh-Taylor',
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
	--initState = 'self-gravitation test 1',
	--initState = 'self-gravitation test 1 spinning',
	--initState = 'self-gravitation test 2',
	--initState = 'self-gravitation test 2 orbiting',
	--initState = 'self-gravitation test 4',
	--initState = 'self-gravitation soup',
	
	-- those designed for SRHD / GRHD:
	--initState = 'relativistic shock reflection',			-- not working.  these initial conditions are constant =P
	--initState = 'relativistic blast wave test problem 1',
	--initState = 'relativistic blast wave test problem 2',
	--initState = 'relativistic blast wave interaction',		-- in 2D this only works with no limiter / lots of dissipation 

	-- states for ideal MHD or two-fluid (not two-fluid-separate)
	--initState = 'Brio-Wu',
	--initState = 'Orszag-Tang',
	--initState = 'MHD rotor',
	--initState = 'spinning magnetic fluid',
	--initState = 'magnetic fluid',
	--initState = '2017 Degris et al',
	--initState = 'that one mhd simulation from youtube',
	--initState = 'GEM challenge', eqnArgs = {useEulerInitState=false},
	
	-- 2002 Dedner
	--initState = '2002 Dedner peak Bx',
	--initState = '2002 Dedner 1D Riemann',
	--initState = '2002 Dedner Shock Reflection',
	--initState = '2002 Dedner 2D Riemann problem',
	--initState = '2002 Dedner Kelvin-Helmholtz',

	-- Maxwell:
	--initState = 'Maxwell default',
	--initState = 'Maxwell scattering around cylinder',
	initState = 'Maxwell scattering around Koch snowflake',
	--initState = 'Maxwell wire',
	
	-- hmm, I think I need a fluid solver for this, not just a Maxwell solver ...
	--initState = 'Maxwell Lichtenberg',	

	-- Maxwell+HD
	--initState = 'two-fluid emhd modified Brio-Wu',
	--initState = 'two-fluid EMHD soliton ion',
	--initState = 'two-fluid EMHD soliton electron',
	--initState = 'two-fluid EMHD soliton maxwell',


	-- Einstein
	--initState = 'Minkowski',
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

	
	--initState = 'black hole - Schwarzschild pseudocartesian',
	
	
	--initState = 'black hole - isotropic',	-- this one has momentum and rotation and almost done with multiple sources.  TODO parameterize

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
	--initState = 'testbed - gauge wave',	-- not working with fd solvers
	--initState = 'testbed - gauge wave - diagonal',
	--initState = 'testbed - linear wave',
	--initState = 'testbed - linear wave - diagonal',
	--initState = 'testbed - Gowdy',


	-- NLS
	--initState = 'Gaussian',
	--initState = 'Ring',
	--initState = 'Oscillatory',
}

-- HD
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='euler'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='euler'})))
--self.solvers:insert(require 'solver.euler-hllc'(args))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='euler'})))

-- still haven't added source terms to this
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='navstokes-wilcox'})))

-- HD - Burgers
-- f.e. and b.e. are working, but none of the r.k. integrators 
-- PLM isn't implemented yet
-- neither is source term / poisson stuff
--self.solvers:insert(require 'solver.euler-burgers'(args))

-- SRHD.  
-- rel blast wave 1 & 2 works in 1D at 256 with superbee flux lim
-- rel blast wave interaction works with superbee flux lim in 1D works at 256, fails at 1024 with float (works with double)
-- 	256x256 double fails with F.E., RK2-Heun, RK2-Ralston, RK2-TVD, RK3, RK4-3/8ths,
-- rel blast wave 1 doesn't work in 64x64. with superbee flux lim
-- rel blast wave 2 with superbee flux lim, Roe solver, works at 64x64 with forward euler
-- 	at 256x256 fails with F.E, RK2, RK2-non-TVD., RK3-TVD, RK4, RK4-TVD, RK4-non-TVD 
--    but works with RK2-Heun, RK2-Ralston, RK2-TVD, RK3, RK4-3/8ths
-- Kelvin-Helmholtz works for all borderes freeflow, float precision, 256x256, superbee flux limiter
--self.solvers:insert(require 'solver.srhd-roe'(args))
--self.solvers:insert(require 'solver.srhd-hll'(args))		-- TODO finishme.  the last piece is 'eigen_forInterface'

-- GRHD
-- this is the solver with plug-ins for ADM metric, 
-- yet doesn't come coupled with any other solver, so it will act just like a srhd solver
--self.solvers:insert(require 'solver.grhd-roe'(args))

-- GRHD+GR
-- here's the GRHD solver with the BSSNOK plugged into it
-- TODO inital condition root finding to make sure the EFE is satisifed
--self.solvers:insert(require 'solver.gr-hd-separate'(args))

-- MHD. 
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

--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='mhd'})))
-- TODO FIXME
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='glm-mhd'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='glm-mhd'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='glm-mhd'})))

-- Maxwell
-- when the state is nonzero, at certain sizes there appear errors in the corners
self.solvers:insert(require 'solver.roe'(table(args, {eqn='maxwell'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='maxwell'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='maxwell'})))

-- GLM Maxwell
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='glm-maxwell'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='glm-maxwell'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='glm-maxwell'})))

-- Maxwell+HD two-fluid electron/ion solver
-- TODO FIXME
-- TODO, with the separate solver, use hll, so the ion, electron, and maxwell all use hll separately
--self.solvers:insert(require 'solver.twofluid-emhd-separate-roe'(args))

-- ...so to try and get around that, here the two are combined into one solver:
-- TODO still needs PLM support
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='twofluid-emhd'})))

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


-- GR

--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm1d_v1'})))
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm1d_v2'})))
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d', eqnArgs={useShift=false}})))
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d', eqnArgs={useShift='MinimalDistortionElliptic'}})))	-- TODO finish me
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d', eqnArgs={useShift='2005 Bona / 2008 Yano'}})))	-- TODO finish me
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d', eqnArgs={useShift='HarmonicShiftCondition-FiniteDifference'}})))	-- breaks, even with b.e. integrator
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='adm3d', eqnArgs={useShift='LagrangianCoordinates'}})))	-- TODO finish me
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='z4'}))) -- TODO finish me

--self.solvers:insert(require 'solver.hll'(table(args, {eqn='adm1d_v1'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='adm1d_v2'})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='adm3d', eqnArgs={useShift=false}})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='adm3d', eqnArgs={useShift='HarmonicShiftCondition-FiniteDifference'}})))
--self.solvers:insert(require 'solver.hll'(table(args, {eqn='z4'}))) -- TODO finish me

--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='adm1d_v1', integrator='backward Euler'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='adm1d_v2', integrator='backward Euler'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='adm3d', integrator='backward Euler'})))
--self.solvers:insert(require 'solver.fdsolver'(table(args, {eqn='z4', integrator='backward Euler'}))) -- TODO finish me


-- the BSSNOK solver works similar to the adm3d for the warp bubble simulation
--  but something gets caught up in the freeflow boundary conditions, and it explodes
-- so I have set constant Minkowski boundary conditions?
-- the BSSNOK solver sometimes explodes / gets errors / nonzero Hamiltonian constraint for forward euler
-- however they tend to not explode with backward euler ... though these numerical perturbations still appear, but at least they don't explode
--self.solvers:insert(require 'solver.bssnok-fd'(args))

-- Z4c finite difference, combining BSSNOK and Z4
--self.solvers:insert(require 'solver.z4c-fd'(args))


--self.solvers:insert(require 'solver.nls'(args))





-- the start of unstructured meshes
--self.solvers:insert(require 'solver.meshsolver'(table(args, {eqn='euler', meshfile='n0012_113-33'})))
-- temp here -- to make sure ordinary solvers still run
--self.solvers:insert(require 'solver.roe'(table(args, {eqn='euler', initState='Sod'})))

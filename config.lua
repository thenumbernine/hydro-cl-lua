--[[
TODO one config per experiment (initial condition + config)
and no more setting config values (boundary, etc) in the init cond file
--]]

local dim = cmdline.dim or 3
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
	fixedDT = cmdline.fixedDT,
	cfl = cmdline.cfl or .6/dim,	-- 1/dim,
	
	--fluxLimiter = cmdline.fluxLimiter or 'superbee',
	--fluxLimiter = 'monotized central',
	fluxLimiter = 'donor cell',
	
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
	
	-- [[ Cartesian
	coord = 'cartesian',
	coordArgs = {vectorComponent='holonomic'},		-- use the coordinate derivatives to represent our vector components (though they may not be normalized)
	--coordArgs = {vectorComponent='anholonomic'},		-- use orthonormal basis to represent our vector components
	--coordArgs = {vectorComponent='cartesian'},			-- use cartesian vector components 
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
				{256,256,1},
				{16,16,16},
			},
			['Intel(R) OpenCL HD Graphics/Intel(R) Gen9 HD Graphics NEO'] = {
				{128,1,1},
				{128,128,1},
				
				-- for 11th WENO (2010 Shen Zha) once we reduce size below 6,6 it breaks
				-- so TODO something about boundary conditions on WENO or something ... maybe an error
				-- other than weno, this works fine with finite volume codes
				--{64,1,1},
				{10,10,10},
			},
		})[platAndDevicesNames]
		-- default size options
		or {
			{256,1,1},
			{256,256,1},
			{32,32,32},
		}
	)[dim],
	boundary = type(cmdline.boundary) == 'table' and cmdline.boundary or {
		xmin = cmdline.boundary or 'freeflow',
		xmax = cmdline.boundary or 'freeflow',
		ymin = cmdline.boundary or 'freeflow',
		ymax = cmdline.boundary or 'freeflow',
		zmin = cmdline.boundary or 'freeflow',
		zmax = cmdline.boundary or 'freeflow',
	},
	--]]
	--[[ cylinder
	coord = 'cylinder',
	--coordArgs = {vectorComponent='holonomic'},		-- use the coordinate derivatives to represent our vector components (though they may not be normalized)
	--coordArgs = {vectorComponent='anholonomic'},		-- use orthonormal basis to represent our vector components
	coordArgs = {vectorComponent='cartesian'},			-- use cartesian vector components 
	mins = cmdline.mins or {0, 0, -1},
	maxs = cmdline.maxs or {1, 2*math.pi, 1},			-- TODO bake the 2π into the coordinate chart so this matches grid/cylinder.  Technically θ→2πθ means it isn't the standard θ variable.  I did this for UI convenience with CFDMesh.
	gridSize = ({
		{128, 1, 1},	-- 1D
		{32, 32, 1},	-- 2D
		{32, 32, 32},	-- 3D
	})[dim],
	boundary = type(cmdline.boundary) == 'table' and cmdline.boundary or {
		-- r
		-- notice, this boundary is designed with cylindrical components in mind, so it will fail with vectorComponent==cartesian 
		xmin=cmdline.boundary or 'freeflow',		
		--xmin=cmdline.boundary or 'cylinderRMin',	-- use this when rmin=0
		--xmin=cmdline.boundary or 'mirror',
		--xmin=cmdline.boundary or {name='mirror', args={restitution=0}},
		xmax=cmdline.boundary or 'freeflow',
		--xmax=cmdline.boundary or 'mirror',
		--xmax=cmdline.boundary or {name='mirror', args={restitution=0}},
		
		-- θ
		ymin=cmdline.boundary or 'periodic',
		ymax=cmdline.boundary or 'periodic',
		
		-- z
		--zmin=cmdline.boundary or {name='mirror', args={restitution=0}},
		--zmax=cmdline.boundary or {name='mirror', args={restitution=0}},
		zmin=cmdline.boundary or 'freeflow',
		zmax=cmdline.boundary or 'freeflow',
	},
	--]]
	--[[ Sphere: r, θ, φ 
	coord = 'sphere',
	--coordArgs = {volumeDim = 3},	-- use higher dimension volume, even if the grid is only 1D to 3D
	mins = cmdline.mins or {0, 0, -math.pi},
	maxs = cmdline.maxs or {8, math.pi, math.pi},
	gridSize = ({
		{160, 1, 1}, -- 1D
		{32, 32, 1}, -- 2D
		{16, 16, 16}, -- 3D
	})[dim],
	boundary = type(cmdline.boundary) == 'table' and cmdline.boundary or {
		xmin=cmdline.boundary or 'sphereRMin',
		xmax=cmdline.boundary or 'freeflow',	--'fixed',
		ymin=cmdline.boundary or 'sphereTheta',
		ymax=cmdline.boundary or 'sphereTheta',
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
	boundary = type(cmdline.boundary) == 'table' and cmdline.boundary or {
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
	--initStateArgs = {v={1e-1,1e-1}},
	
	--initState = 'linear',
	--initState = 'gaussian',
	--initState = 'advect wave',
	--initState = 'sphere',
	initState = 'spiral',
	--initState = 'rarefaction wave',
	--initState = 'Bessel',
	--initState = 'cyclone',
	
	--initState = 'Sod',
	--initState = 'Sod with physical units',
	--initStateArgs = {dim=cmdline.displayDim},
	
	--initState = 'rectangle',
	--initState = 'Sedov',
	--initState = 'Noh',
	--initState = 'implosion',
	--initState = 'Kelvin-Helmholtz',
	--initState = 'Rayleigh-Taylor',	--FIXME ... get initial / static hydro potential working
	--initState = 'Colella-Woodward',
	--initState = 'double mach reflection',
	--initState = 'square cavity',
	--initState = 'shock bubble interaction',		-- with usePLM only works with prim or with athena
	--initState = 'Richmyer-Meshkov',
	--initState = 'radial gaussian',

	-- 2002 Kurganov, Tadmor, "Solution of Two-Dimensional Riemann Problems for Gas Dynamics without Riemann Problem Solvers"
	--initState = 'configuration 1',
	--initState = 'configuration 2',
	--initState = 'configuration 3',
	--initState = 'configuration 4',
	--initState = 'configuration 5',
	--initState = 'configuration 6',
	
	-- states for ideal MHD or two-fluid (not two-fluid-separate)
	--initState = 'Brio-Wu',
	--initState = 'Orszag-Tang',
	--initState = 'MHD rotor',
	--initState = 'GEM challenge', eqnArgs = {useEulerInitState=false},
	--initState = 'spinning magnetic fluid',
	--initState = 'magnetic fluid',
	--initState = '2017 Degris et al',
	--initState = 'that one mhd simulation from youtube',
	--initState = 'spiral with flipped B field',
	
	-- 2002 Dedner
	--initState = '2002 Dedner peak Bx',
	--initState = '2002 Dedner 1D Riemann',
	--initState = '2002 Dedner Shock Reflection',
	--initState = '2002 Dedner 2D Riemann problem',
	--initState = '2002 Dedner Kelvin-Helmholtz',
	
	-- Mara initital conditions
	--initState = 'Mara IsentropicPulse',
	--initState = 'Mara Explosion',
	--initState = 'Mara KelvinHelmholtz',
	--initState = 'Mara SmoothKelvinHelmholtz',
	--initState = 'Mara Shocktube1',
	--initState = 'Mara Shocktube2',
	--initState = 'Mara Shocktube3',
	--initState = 'Mara Shocktube4',
	--initState = 'Mara Shocktube5',
	--initState = 'Mara ContactWave',
	--initState = 'Mara RMHDShocktube1',
	--initState = 'Mara RMHDShocktube2',
	--initState = 'Mara RMHDShocktube3',
	--initState = 'Mara RMHDShocktube4',
	--initState = 'Mara RMHDContactWave',
	--initState = 'Mara RMHDRotationalWave',


	-- self-gravitation tests:
	--initState = 'self-gravitation - Earth',	-- validating units along with self-gravitation.
	--initState = 'self-gravitation test 1',
	--initState = 'self-gravitation test 1 spinning',
	--initState = 'self-gravitation test 2',		--FIXME
	--initState = 'self-gravitation test 2 orbiting',
	--initState = 'self-gravitation test 4',
	--initState = 'self-gravitation soup',	--FIXME


	-- those designed for SRHD / GRHD from Marti & Muller 1998:
	--initState = 'relativistic shock reflection',			-- FIXME.  these initial conditions are constant =P
	--initState = 'relativistic blast wave test problem 1',
	--initState = 'relativistic blast wave test problem 2',
	--initState = 'relativistic blast wave interaction',		-- in 2D this only works with no limiter / lots of dissipation 



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
	--initState = 'Wave-FD Gaussian',
	--initState = 'Wave-FD Bessel',

	-- multi-devices
	multiSlices = {3, 1, 1},
}


if cmdline.solver then self.solvers:insert(require('hydro.solver.'..cmdline.solver)(table(args, cmdline))) return end


-- wave equation


--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='wave'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='wave'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='wave', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='wave', wenoMethod='2008 Borges', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='wave', wenoMethod='2010 Shen Zha', order=5})))
--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='wave'})))

-- wave equation with background spacetime metric
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='wave', eqnArgs={beta={'-y / (r * r)','x / (r * r)','0'}}, wenoMethod='1996 Jiang Shu', order=5})))

--[[ Acoustic black hole 
args.eqnArgs = args.eqnArgs or {}
-- for 2012 Visser A.B. metric, 
args.eqnArgs.alpha = '1' 
-- v0 = e_rHat * A/r + e_phiHat * B/r
-- event horizon = |A|/c
-- ergoregion = sqrt(A^2 + B^2)/c
-- v = 10 (-y,x) / r
--args.eqnArgs.beta = {'-10.*y/r', '10.*x/r'}
-- TODO ... whose coordinates is this in?
--args.eqnArgs.beta = {'-1', '0', '0'}
-- K = -v^i_,i / alpha = 0
-- and it assume gamma_ij == the grid metric gamma_ij
self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='wave'})))
--]]


-- shallow water equations


--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='shallow-water'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='shallow-water'})))


-- compressible Euler equations


--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler'})))

--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='euler', hllCalcWaveMethod='Davis direct bounded'})))	-- this is the default hllCalcWaveMethod
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='euler', hllCalcWaveMethod='Davis direct'})))

--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='euler'})))

--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='euler-hllc', eqn='euler', fluxArgs={hllcMethod=0}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='euler-hllc', eqn='euler', fluxArgs={hllcMethod=1}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='euler-hllc', eqn='euler', fluxArgs={hllcMethod=2}})))

-- NOTICE, these are very accurate with RK4, etc., but incur oscillations with Forward-Euler
-- TODO weno doesn't seem to work with self-gravitation
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2008 Borges', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=5})))

--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=7})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2008 Borges', order=7})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=7})))

-- order=9 and above 2D Sod test, things tend to explode ... maybe I should take away the betaCoeffs denominator (since they get normalized anyways?)
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=9})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2008 Borges', order=9})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=9})))

--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=11})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2008 Borges', order=11})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=11})))

-- hmm, 2D Sod 64x64 RK4 fails at just past 1 second ... 
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=13})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2008 Borges', order=13})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=13})))

-- blows up.  maybe I need an implicit RK scheme...
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=13, integrator='backward Euler'})))

--[[ testing different flux methods
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=5, fluxMethod='Lax-Friedrichs'})))
self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=5, fluxMethod='Roe'})))
-- FIXME: Marquina is broken:
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=5, fluxMethod='Marquina'})))
--]]

-- incompressible...
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler', eqnArgs={incompressible=true}})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', eqnArgs={incompressible=true}, wenoMethod='2010 Shen Zha', order=5})))


-- TODO FIXME 
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='navstokes-wilcox'})))


-- compressible Euler equations - Burgers solver


-- f.e. and b.e. are working, but none of the r.k. integrators
-- PLM isn't implemented yet
-- neither is source term / poisson stuff
--self.solvers:insert(require 'hydro.solver.euler-burgers'(args))


-- special relativistic compressible hydrodynamics


-- TODO FIXME
-- rel blast wave 1 & 2 works in 1D at 256 with superbee flux lim
-- rel blast wave interaction works with superbee flux lim in 1D works at 256, fails at 1024 with float (works with double)
-- 	256x256 double fails with F.E., RK2-Heun, RK2-Ralston, RK2-TVD, RK3, RK4-3/8ths,
-- rel blast wave 1 doesn't work in 64x64. with superbee flux lim
-- rel blast wave 2 with superbee flux lim, Roe solver, works at 64x64 with forward euler
-- 	at 256x256 fails with F.E, RK2, RK2-non-TVD., RK3-TVD, RK4, RK4-TVD, RK4-non-TVD 
--    but works with RK2-Heun, RK2-Ralston, RK2-TVD, RK3, RK4-3/8ths
-- Kelvin-Helmholtz works for all borderes freeflow, float precision, 256x256, superbee flux limiter
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='srhd'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='srhd'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='srhd', wenoMethod='2010 Shen Zha', order=5})))

-- srhd incompressible
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='srhd', eqnArgs={incompressible=true}})))


-- general relativistic compressible hydrodynamics


-- TODO remove calcEigenBasis from hydro/eqn/grhd.cl
-- this is the solver with plug-ins for ADM metric, 
-- yet doesn't come coupled with any other solver, so it will act just like a srhd solver
--self.solvers:insert(require 'hydro.solver.grhd-roe'(args))


-- ideal magnetohydrodynamics


-- with superbee flux lim:  
-- Brio-Wu works in 1D at 256, works in 2D at 64x64 in a 1D profile in the x and y directions.
-- Orszag-Tang with forward Euler integrator fails at 64x64 around .7 or .8
-- 		but works with 'Runge-Kutta 4, TVD' integrator at 64x64
-- 		RK4-TVD fails at 256x256 at just after t=.5
--		and works fine with backwards Euler 
-- when run alongside HD Roe solver, curves don't match (different heat capacity ratios?)
--		but that could be because of issues with simultaneous solvers.
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='mhd'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='mhd'})))

-- explodes with Orszag-Tang
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='mhd', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='mhd', wenoMethod='2010 Shen Zha', order=5})))

--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='mhd'})))
-- eqn.useFixedCh == false is failing
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='glm-mhd'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='glm-mhd'})))
--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='glm-mhd'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='glm-mhd', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='glm-mhd', wenoMethod='2010 Shen Zha', order=5})))


-- Maxwell


--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='maxwell'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='maxwell'})))
--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='maxwell'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='maxwell', wenoMethod='1996 Jiang Shu', order=5})))


-- GLM Maxwell


--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='glm-maxwell'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='glm-maxwell'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='glm-maxwell', wenoMethod='2010 Shen Zha', order=7})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='glm-maxwell', wenoMethod='2010 Shen Zha', order=13})))
--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='glm-maxwell'})))
-- doesn't seem to be accurate with any cfls larger than the explicit integrators...
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='glm-maxwell', integrator='backward Euler', integratorArgs={epsilon=1e-19, verbose=true}, fixedDT=.125})))


-- Maxwell+HD two-fluid electron/ion solver
-- TODO FIXME
-- TODO, with the separate solver, use hll, so the ion, electron, and maxwell all use hll separately
-- TODO I made it even more difficult to implement with the addition of these real_ and cplx_ macros...
--self.solvers:insert(require 'hydro.solver.twofluid-emhd-separate-roe'(args))

-- ...so to try and get around that, here the two are combined into one solver:
-- it is stable, however since all variables are tied together, it integrates them together, which means everything is explicit-integration updated
-- ...which means, with the Maxwell equations waves propagating at the speed of light, that it goes very slow
-- TODO: I suppose I could make this work with my integrator by (1) removing the maxwell terms from the integration variable list and (2) providing a separate operator that updates them implicitly
-- TODO still needs PLM support
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='twofluid-emhd'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='twofluid-emhd', wenoMethod='1996 Jiang Shu', order=9})))	-- exploded...
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='twofluid-emhd', wenoMethod='2010 Shen Zha', order=5})))


-- here's another one: two-fluid emhd with de Donder gauge linearized general relativity
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='twofluid-emhd-lingr'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='twofluid-emhd-lingr', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='twofluid-emhd-lingr', wenoMethod='2010 Shen Zha', order=7})))


-- general relativity


--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='adm1d_v1'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='adm1d_v2'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='adm3d'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='adm3d', eqnArgs={noZeroRowsInFlux=false}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='adm3d', eqnArgs={useShift='MinimalDistortionElliptic'}})))	-- TODO finish me
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='adm3d', eqnArgs={useShift='MinimalDistortionEllipticEvolve'}})))	-- TODO finish me
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='adm3d', eqnArgs={useShift='2005 Bona / 2008 Yano'}})))	-- TODO finish me
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='adm3d', eqnArgs={useShift='HarmonicShiftCondition-FiniteDifference'}})))	-- breaks, even with b.e. integrator
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='adm3d', eqnArgs={useShift='LagrangianCoordinates'}})))	-- TODO finish me
-- TODO bring these up to date with the new normalInfo_t system:
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='z4_2008yano'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='z4'})))

--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='adm1d_v1'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='adm1d_v2'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='adm3d'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='adm3d', eqnArgs={noZeroRowsInFlux=false}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='adm3d', eqnArgs={useShift='HarmonicShiftCondition-FiniteDifference'}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='z4'})))

--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='adm3d', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='adm3d', wenoMethod='2010 Shen Zha', order=7})))

--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='adm1d_v1', integrator='backward Euler'})))
--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='adm1d_v2', integrator='backward Euler'})))
--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='adm3d', integrator='backward Euler'})))
--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='z4_2008yano', integrator='backward Euler'})))



--[[
bssnok is working in 1D-3D Cartesian for RK4 
diverging for non-Cartesian
for spherical, Minkowski init cond, dim=1, gridSize=256, no dissipation, numGhost=2 <=> derivOrder=4, cfl=.4
for RK4 integrator, 1D, range [0, 8] 
	runs indefinitely
for 2D, [40,40] grid, hyperbolic gamma driver shift,
	using the 2017 Ruchlin cfl condition 
		ran until |H|>1 at 0.18491652961003 	
	using the 2008 Alcubierre Hyperbolic Formalism chapter's equations for speed-of-light (alpha sqrt(gamma^ii)):
	using cfl=.2
		runs until |H|>1 at t=3.0379149106647
--]]
--self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn='bssnok-fd-num'})))	-- default shift is HyperbolicGammaDriver
--self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn='bssnok-fd-num', eqnArgs={useShift='none'}})))
--self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn='bssnok-fd-num', eqnArgs={useShift='GammaDriver'}})))

--[[
BSSNOK but with my symbolic CAS generating the math
Generation is really slow and not yet cached.
In spherical on my laptop this is ~1min to do the differentiatin and simpliciations, then ~5min to compile (as opposed to the ~20sec to compile the bssnok-fd-num version).
This one, Minkowski spherical, with no shift is completely smooth.  no oscillations, runs indefinitely.
	always keeping a |H| ~ 1e-10 for vacuum spacetimes
With hyperbolic gamma driver shift it has trouble.
--]]
--self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn='bssnok-fd-sym'})))

-- Z4c finite difference, combining BSSNOK and Z4
-- FIXME something is asymmetric.  watch Theta.  Run warp bubble.
--self.solvers:insert(require 'hydro.solver.z4c-fd'(args))



-- GRHD+GR
-- here's the GRHD solver with the BSSNOK plugged into it
-- TODO inital condition root finding to make sure the EFE is satisifed
--self.solvers:insert(require 'hydro.solver.gr-hd-separate'(args))


-- GR+Maxwell.  params go to the Maxwell solver.

-- This is a Maxwell solver with the extra terms to account for GR.
-- It only simulates Maxwell -- it needs an additional GR solver to fully work correctly.
-- (Without a GR solver it will just operate in Minkowski spacetime.)
-- TODO this should be the einstein-maxwell eqns solver with flat space plugged into the metric
-- TODO rename to einstein-maxwell-roe
-- TODO it's unfinished
--self.solvers:insert(require 'hydro.solver.gr-maxwell-roe'(args))

-- TODO and this should be the einstein-maxwell with an einstein solver plugged in
-- TODO rename to einstein-maxwell-separate
--self.solvers:insert(require 'hydro.solver.gr-em-separate'(args))


-- nonlinear Schrodinger equation
--self.solvers:insert(require 'hydro.solver.nls'(args))

--[[ wave with background metric of acoustic black hole, fourier transform mode, finite difference
-- TODO rename to something more telling, like 'wave-ab-mode-fd'
self.solvers:insert(require 'hydro.solver.wave-fd'(args))
-- I think b.e. is broke
self.solvers:insert(require 'hydro.solver.wave-fd'(table(args, {integrator='backward Euler'})))
-- referencing this integrator doesn't seem to work
--self.solvers:insert(require 'hydro.solver.wave-fd'(table(args, {integrator='Iterative Crank-Nicolson'})))
-- TODO how about a weno-finite-difference solver?
--]]



-- unstructured meshes
-- CFDMesh runs 50x50 at 100 fps (TODO without display)
-- CFDMesh runs 256x256 at 10.5 fps (TODO without display)
-- TODO Hydro 50x50 with fluxLimiter
-- TODO Hydro 50x50 without fluxLimiter
-- hydro-cl GridSolver with fluxLimiter runs 50x50 at 800 fps
-- hydro-cl GridSolver without fluxLimiter runs 50x50 at 2200 fps
-- hydro-cl GridSolver with fluxLimiter runs 256x256 at 60 fps
-- hydro-cl GridSolver without fluxLimiter runs 256x256 at 155 fps
-- hydro-cl MeshSolver runs 50x50 at 2500 fps
-- hydro-cl MeshSolver runs 256x256 at 70 fps (building the mesh took 4.5 minutes =P)
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2d', triangulate=true, size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='p2dfmt', meshfile='n0012_113-33.p2dfmt'}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2dcbrt', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2dcubed', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='cylinder2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='polar2d', size={64, 64}, mins={0,0}}})))

-- TODO FIXME
self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='cube3d', size={8, 8, 8}}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='cylinder3d', size={16, 16, 16}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='sphere3d', size={16, 16, 16}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='torus3d', size={16, 16, 16}}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', eqnArgs={incompressible=true}, mesh={type='image2d', extrude=1, image='blueprints/blueprint.png'}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='hll', eqn='euler', mesh={type='polar2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='euler-hllc', eqn='euler', mesh={type='polar2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='mhd', mesh={type='polar2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='glm-mhd', mesh={type='polar2d', size={64, 64}}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='wave', mesh={type='polar2d', size={64, 64}, mins={0,0,0}}})))

-- rmin=0 without capmin is working
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='mhd', mesh={type='polar2d', size={64, 64}, mins={0,0,0}}})))
-- rmin=0 with capmin is crashing
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='mhd', mesh={type='polar2d', size={64, 64}, mins={0,0,0}, capmin={1,0,0}}})))

-- NEXT BIG TODO
-- * make flux an option in gridsolver as well
-- * make meshsolver and gridsolver separate options
-- * make grid coordinate chart separate of vector component coordinate chart
-- so in the end, the user can choose: 
-- 1) geometry: mesh (w/coordinate mapping) vs grid
-- 2) flux
-- 3) vector component coordinates


-- multi GPU
-- how about 'composite grid' instead of 'chopped up'?
--self.solvers:insert(require 'hydro.solver.choppedup'(table(args, {flux='roe', eqn='euler', subsolverClass=require 'hydro.solver.fvsolver'})))



-- the start of AMR
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler'})))
--self.solvers:insert(require 'hydro.solver.amr'(require 'hydro.solver.fvsolver')(table(args, {flux='roe', eqn='euler'})))



-- [=[ 2013 Baumgarte et al, section IV A 1 example & 2017 Ruchlin, Etienne
local dim = 3
local args = {
	app = self,
	
	eqn = 'bssnok-fd-num', 
	--eqn = 'bssnok-fd-sym', 
	
	eqnArgs = {
		--useShift = 'none',
		--useScalarField = true,	-- needed for the scalar field init cond below
	
		--cflMethod = '2008 Alcubierre',
		--cflMethod = '2013 Baumgarte et al, eqn 32',
		cflMethod = '2017 Ruchlin et al, eqn 53',
	},
	dim = dim,
	integrator = 'Runge-Kutta 4',
	--integrator = 'backward Euler',
	--integratorArgs = {verbose=true},
	cfl = .5,
	
	-- [[
	--coord = 'sphere',
	coord = 'sphere-log-radial',
	mins = {0, 0, 0},
	maxs = {1, math.pi, 2*math.pi},
	gridSize = cmdline.gridSize or ({
		{128, 1, 1},
		{64, 16, 1},
		
		{32, 2, 2},		-- brill-lindquist head-on merger: -- {400, 64, 2},	--{128, 2, 2},
		--{128,32,2},
		--{400, 64, 2},
	})[dim],
	boundary = {
		xmin='sphereRMin',
		
		-- runs a gauge wave until t=...
		--xmax='fixed',		-- 5.3875
		--xmax='freeflow',	-- diverges near rmin after t=60 or so
		--xmax='linear',	-- 13.1875
		xmax='quadratic',	-- 10.6125
		
		ymin='sphereTheta',
		ymax='sphereTheta',
		zmin='periodic',	-- spherePhi is the same as periodic
		zmax='periodic',
	},
	--]]
	--[[
	coord = 'cartesian',
	mins = {-4,-4,-4},
	maxs = {4,4,4},
	gridSize = ({
		{250, 1, 1},
		{40, 40, 1},
		{16, 16, 16},
	})[dim],

	boundary = {
		xmin='freeflow',
		xmax='freeflow',
		ymin='freeflow',
		ymax='freeflow',
		zmin='freeflow',
		zmax='freeflow',
	},
	--]]
	
	--initState = 'Minkowski',	-- TODO sphere-log-radial 
	
	-- TODO look up Teukolsky Phys Rev 26 745 1982 
	--initState = 'pure gauge wave',
	--initState = 'scalar field',
	
	--initState = 'gaussian perturbation',	-- TODO restore this to the 2008 Alcubeirre and 1998 Alcubierre gauge wave examples
	
	--[[
	--initState = 'black hole - boosted Schwarzschild',
	initState = 'black hole - Schwarzschild isotropic - spherical',
	--initState = 'black hole - Brill Lindquist',
	initStateArgs = {
		center = {0,0,0},
		R = .1,
	},
	--]]
	
	--initState = 'black hole - isotropic',
	
	--initState = 'Alcubierre warp bubble',

	--[[
	initState = 'Alcubierre warp bubble',
	initStateArgs = {
		H = -1e-7,	--H = 1e-7,
		sigma = 10,
	},
	--]]
	
	-- only for bssnok-fd-senr
	--initState = 'SENR sphere-log-radial Minkowski',
	initState = 'SENR sphere-log-radial UIUC',
	--initState = 'SENR sphere-log-radial BrillLindquist',
	--initState = 'SENR sphere-log-radial BoostedSchwarzschild',
	--initState = 'SENR sphere-log-radial StaticTrumpet',
}
--self.solvers:insert(require 'hydro.solver.bssnok-fd-pirk'(table(args, {eqn = 'bssnok-fd-num'})))
--self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn = 'bssnok-fd-num'})))
--self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn = 'bssnok-fd-sym'})))
--self.solvers:insert(require 'hydro.solver.bssnok-fd-senr'(args))
if cmdline.bssnok_fd_num then
	self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn = 'bssnok-fd-num'})))
end
if cmdline.bssnok_fd_senr then
	self.solvers:insert(require 'hydro.solver.bssnok-fd-senr'(args))
end
--]=]

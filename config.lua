--[[
TODO one config per experiment (initial condition + config)
and no more setting config values (boundary, etc) in the init cond file
--]]

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
	--integrator = 'backward Euler, CPU',
	--integratorArgs = {verbose=true},
	
	--fixedDT = .0001,
	fixedDT = cmdline.fixedDT,
	cfl = cmdline.cfl or .3/dim,	-- 1/dim,
	
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
	slopeLimiter = 'minmod',
	--slopeLimiter = 'monotized central',
	--slopeLimiter = 'superbee',

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
				{600,1,1},
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
		-- TODO doesn't work
	--coordArgs = {vectorComponent='holonomic'},		-- use the coordinate derivatives to represent our vector components (though they may not be normalized)
		-- TODO doesn't work for rmin=0
	coordArgs = {vectorComponent='anholonomic'},		-- use orthonormal basis to represent our vector components.
	--coordArgs = {vectorComponent='cartesian'},		-- use cartesian vector components 
	mins = cmdline.mins or {0, 0, -1},
	maxs = cmdline.maxs or {1, 2*math.pi, 1},			-- TODO bake the 2π into the coordinate chart so this matches grid/cylinder.  Technically θ→2πθ means it isn't the standard θ variable.  I did this for UI convenience with CFDMesh.
	gridSize = ({
		{128, 1, 1},	-- 1D
		{64, 64, 1},	-- 2D
		{32, 32, 32},	-- 3D
	})[dim],
	boundary = type(cmdline.boundary) == 'table' and cmdline.boundary or {
		-- r
		-- notice, this boundary is designed with cylindrical components in mind, so it will fail with vectorComponent==cartesian 
		--xmin=cmdline.boundary or 'freeflow',		
		xmin=cmdline.boundary or 'cylinderRMin',	-- use this when rmin=0
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
	--coordArgs = {vectorComponent='holonomic'},
	coordArgs = {vectorComponent='anholonomic'},
	--coordArgs = {vectorComponent='cartesian'},
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

	-- TODO separate initConds for each class of equation
	-- this would cohese better with the combined solvers
	-- i.e. a fluid initCond, a Maxwell init-state, and a GR init-state
	-- ... but that means splitting the MHD init-states across M and HD ...
	-- how about just stacking initConds?
	-- and letting each one assign what values it wants.
	-- still to solve -- how do we specify initConds for combined solvers with multiple sets of the same variables (ion/electron, etc)

	-- no initial state means use the first
	--initCond = cmdline.initCond,
	
	-- Euler / SRHD / MHD initial states:
	
	--initCond = 'constant',
	--initCondArgs = {v={1e-1,1e-1}},
	
	--initCond = 'random',
	--initCond = 'linear',
	--initCond = 'gaussian',
	--initCond = 'advect wave',
	--initCond = 'sphere',
	--initCond = 'spiral',
	--initCond = 'rarefaction wave',
	--initCond = 'Bessel',
	--initCond = 'cyclone',
	
	--initCond = 'Sod',
	--initCond = 'Sod with physical units',
	--initCondArgs = {dim=cmdline.displayDim},
	
	--initCond = 'rectangle',
	--initCond = 'Sedov',
	--initCond = 'Noh',
	--initCond = 'implosion',
	--initCond = 'Kelvin-Helmholtz',
	--initCond = 'Rayleigh-Taylor',	--FIXME ... get initial / static hydro potential working
	--initCond = 'Colella-Woodward',
	--initCond = 'double mach reflection',
	--initCond = 'square cavity',
	--initCond = 'shock bubble interaction',		-- with usePLM only works with prim or with athena
	--initCond = 'Richmyer-Meshkov',
	--initCond = 'radial gaussian',

	-- 2002 Kurganov, Tadmor, "Solution of Two-Dimensional Riemann Problems for Gas Dynamics without Riemann Problem Solvers"
	--initCond = 'configuration 1',
	--initCond = 'configuration 2',
	--initCond = 'configuration 3',
	--initCond = 'configuration 4',
	--initCond = 'configuration 5',
	--initCond = 'configuration 6',
	
	-- states for ideal MHD or two-fluid (not two-fluid-separate)
	--initCond = 'Brio-Wu',
	--initCond = 'Orszag-Tang',
	--initCond = 'MHD rotor',
	--initCond = 'GEM challenge', eqnArgs = {useEulerInitState=false},
	--initCond = 'spinning magnetic fluid',
	--initCond = 'magnetic fluid',
	--initCond = '2017 Degris et al',
	--initCond = 'that one mhd simulation from youtube',
	--initCond = 'spiral with flipped B field',
	
	-- 2002 Dedner
	--initCond = '2002 Dedner peak Bx',
	--initCond = '2002 Dedner 1D Riemann',
	--initCond = '2002 Dedner Shock Reflection',
	--initCond = '2002 Dedner 2D Riemann problem',
	--initCond = '2002 Dedner Kelvin-Helmholtz',
	
	-- Mara initital conditions
	--initCond = 'Mara IsentropicPulse',
	--initCond = 'Mara Explosion',
	--initCond = 'Mara KelvinHelmholtz',
	--initCond = 'Mara SmoothKelvinHelmholtz',
	--initCond = 'Mara Shocktube1',
	--initCond = 'Mara Shocktube2',
	--initCond = 'Mara Shocktube3',
	--initCond = 'Mara Shocktube4',
	--initCond = 'Mara Shocktube5',
	--initCond = 'Mara ContactWave',
	--initCond = 'Mara RMHDShocktube1',
	--initCond = 'Mara RMHDShocktube2',
	--initCond = 'Mara RMHDShocktube3',
	--initCond = 'Mara RMHDShocktube4',
	--initCond = 'Mara RMHDContactWave',
	--initCond = 'Mara RMHDRotationalWave',


	-- self-gravitation tests:
	--initCond = 'self-gravitation - Earth',	-- validating units along with self-gravitation.
	--initCond = 'self-gravitation test 1',
	--initCond = 'self-gravitation test 1 spinning',
	--initCond = 'self-gravitation test 2',		--FIXME
	--initCond = 'self-gravitation test 2 orbiting',
	--initCond = 'self-gravitation test 4',
	--initCond = 'self-gravitation soup',	--FIXME


	-- those designed for SRHD / GRHD from Marti & Muller 1998:
	--initCond = 'relativistic shock reflection',			-- FIXME.  these initial conditions are constant =P
	--initCond = 'relativistic blast wave test problem 1',
	--initCond = 'relativistic blast wave test problem 2',
	--initCond = 'relativistic blast wave interaction',		-- in 2D this only works with no limiter / lots of dissipation 



	-- Maxwell:
	--initCond = 'Maxwell default',
	--initCond = 'Maxwell empty waves',
	--initCond = 'Maxwell scattering around cylinder',
	--initCond = 'Maxwell scattering around pyramid',
	--initCond = 'Maxwell scattering around square',
	--initCond = 'Maxwell scattering around Koch snowflake',
	--initCond = 'Maxwell wire',
	--initCond = 'Maxwell transverse waves',
	--initCond = 'Maxwell charged particle',

	-- hmm, I think I need a fluid solver for this, not just a Maxwell solver ...
	--initCond = 'Maxwell Lichtenberg',	

	-- Maxwell+HD
	--initCond = 'two-fluid emhd modified Brio-Wu',
	--initCond = 'two-fluid EMHD soliton ion',
	--initCond = 'two-fluid EMHD soliton electron',
	--initCond = 'two-fluid EMHD soliton maxwell',


	-- Einstein
	--initCond = 'Minkowski',
	--initCond = 'gaussian perturbation',
	--initCond = 'plane gauge wave',


	--initCond = 'Alcubierre warp bubble',
	
	--initCondArgs = {R=.5, sigma=8, speed=.1},	-- sub-luminal
	
	--initCondArgs = {R=.5, sigma=8, speed=1.1},		-- super-luminal 1.1x
	-- ... works with
	--	size=64x64 solver=adm3d int=fe plm=athena ctu
	--  size=64x64 solver=adm3d int=fe plm=athena
	--  size=64x64 solver=adm3d int=fe flux-limiter=superbee
	
	--initCondArgs = {R=.5, sigma=8, speed=2},		-- super-luminal 2x
	-- ... works with
	--  size=64x64 solver=adm3d int=fe flux-limiter=superbee
	
	--initCondArgs = {R=.5, sigma=8, speed=10},		-- super-luminal 10x
	--  size=64x64 solver=roe eqn=adm3d int=fe flux-limiter=superbee ... eventually explodes
	--  size=64x64 solver=roe eqn=bssnok int=be flux-limiter=superbee ... eventually explodes as well
	--  size=128x128 solver=hll eqn=adm3d int=fe flux-limiter=superbee ... runs for a really long time 

	
	--initCond = 'black hole - Schwarzschild',
	
	
	--initCond = 'black hole - isotropic',	-- this one has momentum and rotation and almost done with multiple sources.  TODO parameterize
	
	--initCond = 'black hole - SENR/NumPy',

	--[[ single black hole, stationary
	initCondArgs = {
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
	initCondArgs = {
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
	initCondArgs = {
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
	initCondArgs = {
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
	initCondArgs = {
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
	initCondArgs = {
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
	initCondArgs = {
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
	initCondArgs = {
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


	--initCond = 'stellar model',
	--initCond = '1D black hole - wormhole form',

	
	--initCond = 'Gowdy waves',
	--initCond = 'testbed - robust',	-- not working with fv solvers
	--initCond = 'testbed - gauge wave',	-- not working with forward-euler finite-difference solvers
	--initCond = 'testbed - gauge wave - diagonal',
	--initCond = 'testbed - linear wave',
	--initCond = 'testbed - linear wave - diagonal',
	--initCond = 'testbed - Gowdy',


	-- NLS
	--initCond = 'Gaussian',
	--initCond = 'Ring',
	--initCond = 'Oscillatory',
	---- piggybacking for my wave-finite-difference here: ----
	initCond = 'Wave-FD Gaussian',
	--initCond = 'Wave-FD Bessel',

	-- multi-devices
	multiSlices = {cmdline.multiSlices or 3, 1, 1},
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
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='shallow-water', wenoMethod='1996 Jiang Shu', order=5})))


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

-- Navier-Stokes-Wilcox:
-- TODO FIXME not done yet 
-- needs 1) arbitrary normal calculations, 2) module based codegen
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='navstokes-wilcox'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='navstokes-wilcox', eqnArgs={incompressible=true}})))


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
--self.solvers:insert(require 'hydro.solver.grhd-fvsolver'(table(args, {flux='roe'})))


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

--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='mhd', eqnArgs={incompressible=true}})))


-- eqn.useFixedCh == false is failing
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='glm-mhd'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='glm-mhd'})))
--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='glm-mhd'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='glm-mhd', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='glm-mhd', wenoMethod='2010 Shen Zha', order=5})))

--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='glm-mhd', eqnArgs={incompressible=true}})))


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

-- [[ wave with background metric of acoustic black hole, fourier transform mode, finite difference
-- TODO rename to something more telling, like 'wave-ab-mode-fd'
--self.solvers:insert(require 'hydro.solver.wave-fd'(args))
-- I think b.e. is broke
--self.solvers:insert(require 'hydro.solver.wave-fd'(table(args, {integrator='backward Euler'})))
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
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2d', size={16, 16}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2d', triangulate=true, size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='p2dfmt', meshfile='n0012_113-33.p2dfmt'}})))	-- TODO needs boundary conditions
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2dcbrt', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2dcubed', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='cylinder2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='cube3d', size={16, 16, 16}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='cylinder3d', size={16, 16, 16}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='sphere3d', size={16, 16, 16}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='torus3d', size={16, 16, 16}}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', eqnArgs={incompressible=true}, mesh={type='image2d', extrude=1, image='blueprints/blueprint.png'}})))

-- polar, rmin=0, duplicated vtxs in the middle
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='polar2d', size={8, 8}}})))

-- polar, rmin=0, single vtx in the middle
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='polar2d', size={8, 8}, capmin={1,0}}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='hll', eqn='euler', mesh={type='polar2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='euler-hllc', eqn='euler', mesh={type='polar2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='mhd', mesh={type='polar2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='glm-mhd', mesh={type='polar2d', size={64, 64}}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='wave', mesh={type='polar2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='hll', eqn='wave', mesh={type='polar2d', size={64, 64}}})))

-- rmin=0 without capmin is working
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='mhd', mesh={type='polar2d', size={64, 64}}})))
-- rmin=0 with capmin is crashing
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='mhd', mesh={type='polar2d', size={64, 64}, capmin={1,0}}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='glm-mhd', mesh={type='torus3d', size={16, 16, 16}}})))

-- NEXT BIG TODO
-- * make meshsolver and gridsolver separate options
-- * make grid coordinate chart separate of vector component coordinate chart
-- so in the end, the user can choose: 
-- 1) geometry: mesh (w/coordinate mapping) vs grid
--		The only dif between grid and mesh solver is that grid has easy n'th order finite difference stencils, thanks to its parameterization
-- 		I can unify this with mesh if I store per-cell the d/dx^i(x + j*dx^i) in the i'th face direction to the j'th step.  
--		Then I could rewrite the current higher order PLM/WENO stuff for the MeshSolver
-- 2) flux
-- 3) vector component coordinates


-- multi GPU
-- how about 'composite grid' instead of 'chopped up'?
--self.solvers:insert(require 'hydro.solver.choppedup'(table(args, {flux='roe', eqn='euler', subsolverClass=require 'hydro.solver.fvsolver'})))



-- the start of AMR
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler'})))
--self.solvers:insert(require 'hydro.solver.amr'(require 'hydro.solver.fvsolver')(table(args, {flux='roe', eqn='euler'})))


-- [=[ reproducing 2009 Alic, Bona, Bona-Casas"Towards a gauge-polyvalent numerical relativity code"
local dim = 1
local args = {
	app = self,
	
	--eqn = 'adm3d',
	eqn = 'z4',
	
	--integrator = 'forward Euler',
	integrator = 'Runge-Kutta 3, TVD',	-- p.20, eqn B.1
	dim = dim,
	cfl = .5/dim,			-- no mention of cfl or timestep ...
	fluxLimiter = cmdline.fluxLimiter or 'donor cell',	-- I didn't read much into what kind of flux/slope limiter was used

	--[[
	coord = 'cartesian',
	mins = {-20,-20,-20},
	maxs = {20,20,20},
	gridSize = ({
		{200, 1, 1},
		{40, 40, 1},
		{8, 8, 8},
	})[dim],
	boundary = {
		xmin = 'quadratic',
		xmax = 'quadratic',
		ymin = 'quadratic',
		ymax = 'quadratic',
		zmin = 'quadratic',
		zmax = 'quadratic',
	},
	--]]
	--[[
	coord = 'sphere',
	coordArgs = {
		--vectorComponent = 'holonomic',	-- TODO this is techically the case, but there may be bugs in this.
		vectorComponent = 'anholonomic',
		--vectorComponent = 'cartesian',	-- adm3d isn't designed for cartesian / mesh / arbitrary normals yet
	},
	mins = {0, 0, 0},
	maxs = {
		cmdline.rmax or 20,	-- bottom of p.11: r(max)=20M
		math.pi,
		2*math.pi,
	},
	gridSize = cmdline.gridSize or ({
		{200, 1, 1},		-- bottom of p.11: h = 0.1M
		{64, 16, 1},
		{32, 2, 2},
	})[dim],
	-- boundary described at top of p.21 ... not sure what they mean
	-- advect by maximum propagation speed ... ?
	boundary = {
		xmin='sphereRMin',
		xmax='quadratic',
		ymin='sphereTheta',
		ymax='sphereTheta',
		zmin='periodic',
		zmax='periodic',
	},
	--]]
	-- [[
	coord = 'sphere-sinh-radial',
	coordArgs = {
		-- TODO sort this out
		-- TODO do I have mem write / unwritten vars in "holonomic"?  cuz there seem to be errors that persist past reset()
		-- TODO move cell_area and cell_volume calcs into cell_t fields 
		vectorComponent = 'cartesian',
		--vectorComponent = 'holonomic',	-- our tensor components are holonomic ... except the partial / 1st order state variables, like a_k, d_kij
		--vectorComponent = 'anholonomic',	-- ... these settings also influence the finite volume area/volume calculations (in terms of the vector components) ... 
		
		-- [==[ the paper uses this remapping parameters (eqn 32):
		-- Alic et al: 		R = L sinh(r / L)
		-- comparing Alic's sinh remapping to SENR's remapping:
		-- Etienne et al: 	rho = A sinh(rho / w) / sinh(1/w)
		-- equate R=rho to get:
		-- L sinh(r / L) = A sinh(rho / w) / sinh(1/w)
		-- so with the SENR sinh-remapping function, we get:
		-- w = L, A = L sinh(1/w)
		sinh_w = 1.5,
		amplitude = 1.5 * math.sinh(1 / 1.5),
		--]==]
	},
	mins = {0, 0, 0},
	maxs = {
		20,	-- p.13 says R = 20 M <=> r ~ 463000 M
		-- for L = 1.5, rho(r) = L * sinh(r / L), we find rho(20) = 463078.22018288 ( as p.13 says )
		math.pi,
		2*math.pi,
	},
	gridSize = cmdline.gridSize or ({
		{200, 1, 1},
		{80, 80, 1},
		{16, 8, 8},
	})[dim],
	boundary = {
		xmin='sphereRMin',
		xmax='quadratic',
		ymin='sphereTheta',
		ymax='sphereTheta',
		zmin='periodic',
		zmax='periodic',
	},
	--]]
	--[[ sphere-sinh-radial but with SENR parameters ... for SENR init conds
	coord = 'sphere-sinh-radial',
	coordArgs = {
		vectorComponent = 'cartesian',
		--vectorComponent = 'holonomic',
		--vectorComponent = 'anholonomic',
		sinh_w = .15,
		amplitude = 1000,
	},
	mins = {0, 0, 0},
	maxs = {1, math.pi, 2*math.pi,},
	gridSize = cmdline.gridSize or ({
		{200, 1, 1},
		{80, 80, 1},
		{16, 8, 8},
	})[dim],
	boundary = {
		xmin='sphereRMin',
		xmax='quadratic',
		ymin='sphereTheta',
		ymax='sphereTheta',
		zmin='periodic',
		zmax='periodic',
	},
	--]]

	--initCond = 'Minkowski',
	--initCond = 'gaussian perturbation',
	--initCond = 'plane gauge wave',
	initCond = 'SENR UIUC',
	--initCond = 'SENR BrillLindquist',
	--[[
	-- TODO since converting this to useBSSNVars, it doesn't work for cartesian anymore ...
	initCond = 'black hole - isotropic',	-- this one has momentum and rotation and almost done with multiple sources.  TODO parameterize
	initCondArgs = {
		bodies = {
			R = 2,
			P_u = {0,0,0},
			S_u = {0,0,0},
			pos = {0,0,0},
		}
	},
	--]]

	flux = 'hll',
}
-- comparing hll solvers
if cmdline['2009Alic-adm'] then
	self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {eqn='adm3d'})))
end
if cmdline['2009Alic-z4'] then
	self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {eqn='z4'})))
end
if cmdline['2009Alic-bssnok-fd-senr'] then
	self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn='bssnok-fd-senr'})))
end
if cmdline['2009Alic-bssnok-fd-num'] then
	self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn='bssnok-fd-num'})))
end
if cmdline['2009Alic-bssnok-fd-sym'] then
	self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn='bssnok-fd-sym'})))
end
--]=]




-- [=[ 2013 Baumgarte et al, section IV A 1 example & 2017 Ruchlin, Etienne
local dim = 3
local args = {
	app = self,
	
	eqn = 'bssnok-fd-num', 
	--eqn = 'bssnok-fd-sym', 
	--eqn = 'bssnok-fd-senr', 
	
	eqnArgs = {
		--useShift = 'none',
		--useScalarField = true,	-- needed for the scalar field init cond below
	
		--cflMethod = '2008 Alcubierre',
		--cflMethod = '2013 Baumgarte et al, eqn 32',
		cflMethod = '2017 Ruchlin et al, eqn 53',
	},
	dim = dim,
	
	integrator = cmdline.integrator or 'Runge-Kutta 4',
	--integrator = cmdline.integrator or 'backward Euler',
	--integrator = cmdline.integrator or 'backward Euler, CPU',	-- debugging.   seems that, for grid sizes too small, B.E. GPU fails.  i think because the reduce() gpu function isn't set up for lower bounds of buffer sizes.
	--integratorArgs = {verbose=true},
	cfl = .5,

	--[[
	coord = 'cartesian',
	mins = {-3,-3,-3},
	maxs = {3,3,3},
	gridSize = ({
		{250, 1, 1},
		{40, 40, 1},
		{8, 8, 8},
	})[dim],
	boundary = {
		--xmin='freeflow',
		--xmax='freeflow',
		--ymin='freeflow',
		--ymax='freeflow',
		--zmin='freeflow',
		--zmax='freeflow',
		xmin = 'quadratic',
		xmax = 'quadratic',
		ymin = 'quadratic',
		ymax = 'quadratic',
		zmin = 'quadratic',
		zmax = 'quadratic',
	},
	--]]
	--[[
	coord = 'sphere',
	coordArgs = {
		-- this isn't really used since bssn is a finite-difference solver, so just pick the one that has the least complications.
		-- but if you want to compare this to finite volume then you should use anholonomic or cartesian
		--vectorComponent = 'holonomic',
		vectorComponent = 'anholonomic',
		--vectorComponent = 'cartesian',
	},
	-- mind you, these mins/maxs correlate with SENR
	-- however, this would put the mid phi at pi, which puts the graph on the x- side of the xy plane ... not the x+ side as it would if phi-mid was equal to 0
	mins = {0, 0, 0},
	maxs = {
		-- 2015 Baumgarte et al, spherical coordinates (not sinh-remapped), PIRK uses this rmax:
		cmdline.rmax or 12.8,--24 M1, but I'm using M1 = M2 = 0.5

		math.pi,
		2*math.pi,
	},
	gridSize = cmdline.gridSize or ({
		{128, 1, 1},
		{64, 16, 1},
		
		-- N x 2 x 2:
		{32, 2, 2},
		--{80, 80, 2},
		--{128, 2, 2},
		--{128, 32, 2},
		--{400, 64, 2},
	
		-- 80N x 40N x 2N
		--{160, 80, 4},
	
		-- Brill-Lindquist head-on merger:2017 Ruchlin, Etienne, section 3, 2 paragraphs after eqn 70: 
		--{400, 64, 2},
	
		-- 2015 Baumgarte et al, head-on collision: 128N, 48N, 2
		--{128, 48, 2},
	
		-- SENR PIRK sphere 'agrees with Baumgarte' grid size
		--{64,32,32}, -- seems to be running fine with -num, rk4, sphere, UIUC
		--{16,8,8},
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
	-- [[
	coord = 'sphere-sinh-radial',
	coordArgs = {
		--vectorComponent = 'holonomic',
		vectorComponent = 'anholonomic',
		--vectorComponent = 'cartesian',
		-- SENR uses these parameters:
		amplitude = 1000,
		sinh_w = .15
	},
	mins = {0, 0, 0},
	maxs = {
		1,	-- 2017 Ruchlin et al, rely on coordinate chart to remap to rmax
		math.pi,
		2*math.pi,
	},
	gridSize = cmdline.gridSize or ({
		{128, 1, 1},
		{64, 16, 1},
		
		-- N x 2 x 2:
		--{32, 2, 2},		-- SENR sphere-sinh-radial uses this by default
		--{80, 80, 2},		-- this works well for BrillLindquist sphere-sinh-radial when viewing the xz slice
		--{128, 2, 2},
		--{128, 32, 2},
		--{400, 64, 2},
		
		{200, 2, 2},
	
		-- 80N x 40N x 2N
		--{160, 80, 4},
	
		-- Brill-Lindquist head-on merger:2017 Ruchlin, Etienne, section 3, 2 paragraphs after eqn 70: 
		--{400, 64, 2},
	
		-- 2015 Baumgarte et al, head-on collision: 128N, 48N, 2
		--{128, 48, 2},
	
		--{16,8,8},
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


	--initCond = 'Minkowski',	-- TODO sphere-sinh-radial 
	
	-- TODO look up Teukolsky Phys Rev 26 745 1982 
	--initCond = 'pure gauge wave',
	--initCond = 'scalar field',
	
	--initCond = 'gaussian perturbation',	-- TODO restore this to the 2008 Alcubeirre and 1998 Alcubierre gauge wave examples
	
	--[[
	--initCond = 'black hole - boosted Schwarzschild',
	initCond = 'black hole - Schwarzschild isotropic - spherical',
	--initCond = 'black hole - Brill Lindquist',
	initCondArgs = {
		center = {0,0,0},
		R = .1,
	},
	--]]
	
	--initCond = 'black hole - isotropic',
	
	--initCond = 'Alcubierre warp bubble',

	--[[
	initCond = 'Alcubierre warp bubble',
	initCondArgs = {
		H = -1e-7,	--H = 1e-7,
		sigma = 10,
	},
	--]]
	
	--initCond = 'Minkowski',
	--initCond = 'SENR Minkowski',
	initCond = 'SENR UIUC',
	--initCond = 'SENR BrillLindquist',
	--initCond = 'SENR BoostedSchwarzschild',
	--initCond = 'SENR StaticTrumpet',
}
--self.solvers:insert(require 'hydro.solver.bssnok-fd-pirk'(table(args, {eqn = 'bssnok-fd-num'})))	-- requires extra PIRK kernels to be defined in the eqn file
--self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn = 'bssnok-fd-num'})))
--self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn = 'bssnok-fd-sym'})))
--self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn = 'bssnok-fd-senr'})))
if cmdline['bssnok-fd-num'] then
	self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn = 'bssnok-fd-num'})))
end
if cmdline['bssnok-fd-sym'] then
	self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn = 'bssnok-fd-sym'})))
end
if cmdline['bssnok-fd-senr'] then
	self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn = 'bssnok-fd-senr'})))
end

if cmdline['bssnok-fd-num-pirk'] then
	self.solvers:insert(require 'hydro.solver.bssnok-fd-pirk'(table(args, {eqn = 'bssnok-fd-num'})))
end
if cmdline['bssnok-fd-senr-pirk'] then
	self.solvers:insert(require 'hydro.solver.bssnok-fd-pirk'(table(args, {eqn = 'bssnok-fd-senr'})))
end

if cmdline['adm3d-roe'] then
	self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='adm3d'})))
end
if cmdline['adm3d-hll'] then
	self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='adm3d'})))
end
--]=]

--[[
TODO one config per experiment (initial condition + config)
and no more setting config values (boundary, etc) in the init cond file
--]]
local constants = require 'hydro.constants'
local materials = require 'hydro.materials'

local dim = cmdline.dim or 2
local args = {
	app = self,
	dim = dim,
	eqn = cmdline.eqn,
	flux = cmdline.flux,

	--integrator = cmdline.integrator or 'forward Euler',
	--integrator = 'Iterative Crank-Nicolson',
	--integrator = 'Runge-Kutta 2',
	--integrator = 'Runge-Kutta 2 Heun',
	--integrator = 'Runge-Kutta 2 Ralston',
	--integrator = 'Runge-Kutta 3',
	integrator = 'Runge-Kutta 4',
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

	-- with Kelvin-Helmholts, this will explode even at .5/(dim=2), but runs safe for .3/(dim=2)
	cfl = cmdline.cfl or .3/dim,

	fluxLimiter = cmdline.fluxLimiter or 'superbee',
	--fluxLimiter = 'monotized central',
	--fluxLimiter = 'donor cell',		-- same as turning fluxlimiter off ... you have to turn fluxlimiter off to use plm

	-- piecewise-linear slope limiter
	-- TODO rename this to 'calcLR' or something
	--									-- min div v for gridSize={1024} cfl=.3 Sod mirror at t=0.5:
	--									-- -191 = no plm, superbee flux limiter
	--									-- -184 = no plm, monotized central flux limiter
	--usePLM = 'piecewise-constant',	-- -84		degenerate case.  don't use this, instead just disable usePLM, or else this will allocate more memory / run more functions.
	--usePLM = 'plm-cons',				-- -190
	--usePLM = 'plm-cons-alone',		-- -177
	--usePLM = 'plm-prim-alone',		-- -175
	--usePLM = 'plm-eig',				-- -88		\
	--usePLM = 'plm-eig-prim',			-- -88		 - these have less sharp shock wave in Sod than the non-eig ones
	--usePLM = 'plm-eig-prim-ref',		-- -28 		/
	--usePLM = 'plm-athena',			-- -40		based on Athena.  most accurate from 1D sod tests atm
	--usePLM = 'ppm-wip',				-- 			FIXME one more attempt to figure out all the PLM stuff, based on 2017 Zingale
	--usePLM = 'weno',					-- 			TODO make WENO one of these 'usePLM' methods. rename it to 'construct LR state method' or something.  then use CTU with WENO.  or can we, since even the CTU method should use the re-linear-projection ... i should just have these separate plm methods as separate functions ...

	-- only enabled for certain usePLM methods
	--slopeLimiter = 'minmod',
	--slopeLimiter = 'monotized central',
	--slopeLimiter = 'superbee',

	-- this is functional without usePLM, but doing so falls back on the cell-centered buffer, which with the current useCTU code will update the same cell twice from different threads
	-- TODO this seems to introduce more diagonal waves for SRHD
	--useCTU = true,

	-- [[ Cartesian
	coord = 'cartesian',
	--coordArgs = {vectorComponent='holonomic'},		-- use the coordinate derivatives to represent our vector components (though they may not be normalized)
	--coordArgs = {vectorComponent='anholonomic'},		-- use orthonormal basis to represent our vector components
	coordArgs = {vectorComponent='cartesian'},			-- use cartesian vector components
	mins = cmdline.mins or {-1, -1, -1},
	maxs = cmdline.maxs or {1, 1, 1},

	-- 256^2 = 2^16 = 2 * 32^3
	gridSize = cmdline.gridSize or (
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

			-- latest CL platform/device name on my Intel HD 520 on ubuntu
			-- Intel(R) OpenCL/Intel(R) HD Graphics 520
			-- Intel(R) OpenCL HD Graphics/Intel(R) Gen9 HD Graphics NEO
			-- Intel(R) OpenCL HD Graphics/Intel(R) Gen9 HD Graphics NEO
			-- Intel(R) OpenCL HD Graphics/Intel(R) Graphics Gen9 [0x1916]
			['Intel(R) OpenCL HD Graphics/Intel(R) HD Graphics 520 [0x1916]'] = {
				{4096,1,1},
				{64,64,1},
				{32,32,32},
			},

			-- 5600M with device=gfx902 to work with gl_sharing:
			['AMD Accelerated Parallel Processing/gfx902'] = {
				{4096,1,1},
				{256,256,1},
				{32,32,32},
			},

			-- 5600M with device=gfx1010, which doesn't work with gl_sharing, but runs much faster:
			['AMD Accelerated Parallel Processing/gfx1010'] = {
				{4096,1,1},
				{256,256,1},
				{32,32,32},
			},

			-- 5600M when I'm too lazy to only specify just one device.
			-- it still just uses the first, which is the gfx1010
			-- for 2D:
			--	256x256 runs at about 270 fps
			--	1500x1500 runs at about 20 fps
			-- 4000x4000 uses 4362211620 bytes
			-- from then on, any bigger tends to segfault somewhere after 'randomizing UBuf...'
			['AMD Accelerated Parallel Processing/gfx1010/gfx902'] = {
				{256,1,1},
				{256,256,1},
				{32,32,32},
			},

			-- and on linux ...
			['AMD Accelerated Parallel Processing/gfx1010:xnack-/gfx902:xnack-'] = {
				{256,1,1},
				{256,256,1},
				{48,48,48},
			},
			['CPU debug implementation/CPU debug implementation'] = {
				{256,1,1},
				{64,64,1},
				{32,32,32},
			},
		})[platAndDevicesNames]
		-- default size options
		or {
			{1024,1,1},
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
		-- TODO explodes
	--coordArgs = {vectorComponent='holonomic'},		-- use the coordinate derivatives to represent our vector components (though they may not be normalized)
		-- TODO subtle innacuracy drift
	--coordArgs = {vectorComponent='anholonomic'},		-- use orthonormal basis to represent our vector components.
		-- TODO works but does curvilinear boundaries support cartesian vector components?
		-- for Euler equations / constants / vel = 0.1 e_x + 0.1 e_y the density drifts by 0.3% in the first second
	coordArgs = {vectorComponent='cartesian'},			-- use cartesian vector components
	mins = cmdline.mins or {0, 0, -1},
	maxs = cmdline.maxs or {.5, 2*math.pi, 1},			-- TODO bake the 2π into the coordinate chart so this matches grid/cylinder.  Technically θ→2πθ means it isn't the standard θ variable.  I did this for UI convenience with CFDMesh.
	gridSize = ({
		{128, 1, 1},	-- 1D
		{32, 128, 1},	-- 2D
		{32, 32, 32},	-- 3D
	})[dim],
	boundary = type(cmdline.boundary) == 'table' and cmdline.boundary or {
		-- r
		-- notice, this boundary is designed with cylindrical components in mind, so it will fail with vectorComponent==cartesian
		--xmin=cmdline.boundary or 'freeflow',

		-- solver=fvsolver eqn=wave flux=roe initCond=constant initCondArgs={rho=1, v={.1,.1}}
		-- after 1s: U Pi=[0.98602808771532 1.0004558819272 ±0.0022951109994 1.0155402211526]
		-- the average might look like its further from the init value of rho=1, but it is oscillating around rho=1 and not linearly divering (like xmin='none' causes)
		-- but the error is higher, though it doesn't drift as far
		xmin=cmdline.boundary or 'cylinderRMin',	-- use this when rmin=0

		-- solver=fvsolver eqn=wave flux=roe initCond=constant initCondArgs={rho=1, v={.1,.1}}
		-- after 1s: U Pi=[1.0002535672563 1.0003601365045 ±0.00010068116618478 1.0007179873631]
		-- so the error stays lower, but the average drifts upward
		--xmin=cmdline.boundary or 'none',

		--xmin=cmdline.boundary or 'mirror',
		--xmin=cmdline.boundary or {name='mirror', args={restitution=0}},
		--xmax=cmdline.boundary or 'freeflow',
		xmax=cmdline.boundary or 'freeflow',
		--xmax=cmdline.boundary or {name='mirror', args={restitution=0}},

		-- θ
		ymin=cmdline.boundary or 'periodic',
		ymax=cmdline.boundary or 'periodic',

		-- z
		--zmin=cmdline.boundary or {name='mirror', args={restitution=0}},
		--zmax=cmdline.boundary or {name='mirror', args={restitution=0}},
		zmin=cmdline.boundary or 'freeflow',
		zmax=cmdline.boundary or 'freeflow',
		--zmin=cmdline.boundary or 'mirror',
		--zmax=cmdline.boundary or 'mirror',
	},
	--]]
	--[[ Sphere: r, θ, φ
	coord = 'sphere',
	--coordArgs = {volumeDim = 3},	-- use higher dimension volume, even if the grid is only 1D to 3D
	--coordArgs = {vectorComponent='holonomic'},
	--coordArgs = {vectorComponent='anholonomic'},
	coordArgs = {vectorComponent='cartesian'},
	mins = cmdline.mins or {0, 0, -math.pi},
	maxs = cmdline.maxs or {1, math.pi, math.pi},
	gridSize = ({
		{160, 1, 1}, -- 1D
		{32, 32, 1}, -- 2D
		{32, 32, 32}, -- 3D
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
	coordArgs = {vectorComponent='cartesian'},
	-- hmm, right now sphere's variables change per-dimension used ...
	mins = cmdline.mins or {0, 0, 0},
	maxs = cmdline.maxs or {2*math.pi, 2*math.pi, 1},
	gridSize = cmdline.gridSize or {16,16,16},
	boundary = type(cmdline.boundary) == 'table' and cmdline.boundary or {
		xmin=cmdline.boundary or 'periodic',
		xmax=cmdline.boundary or 'periodic',
		ymin=cmdline.boundary or 'periodic',
		ymax=cmdline.boundary or 'periodic',
		zmin=cmdline.boundary or 'mirror',
		zmax=cmdline.boundary or 'mirror',
	},
	--]]
	--[[ cylinder as toroid
	coord = 'cylinder',
	coordArgs = {vectorComponent='cartesian'},
	--coordArgs = {vectorComponent='anholonomic'},
	mins = cmdline.mins or {.5, 0, -.25},
	maxs = cmdline.maxs or {1, 2*math.pi, .25},			-- TODO bake the 2π into the coordinate chart so this matches grid/cylinder.  Technically θ→2πθ means it isn't the standard θ variable.  I did this for UI convenience with CFDMesh.
	gridSize = ({
		{128, 1, 1},	-- 1D
		{64, 64, 1},	-- 2D
		{32, 32, 32},	-- 3D
	})[dim],
	boundary = type(cmdline.boundary) == 'table' and cmdline.boundary or {
		-- r
		xmin=cmdline.boundary or 'mirror',
		xmax=cmdline.boundary or 'mirror',
		-- θ
		ymin=cmdline.boundary or 'periodic',
		ymax=cmdline.boundary or 'periodic',
		-- z
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
	--initCondArgs = {v={1,0}},
	--initCondArgs = {v={1e-1,1e-1}},

	--initCond = 'random',
	--initCond = 'linear',
	--initCond = 'advect wave',
	--initCond = 'advect gaussian',	-- TODO fix the default case
	--[[ 1D test case
	initCondArgs = {
		rho0 = 1,
		rho1 = 3,
		P0 = 1,
		u0 = .1,
		v0 = 0,
		x0 = -.5,
		y0 = 0,
		z0 = 0,
	},
	--]]
	--[[ 2D test case
	initCondArgs = {
		rho0 = 1,
		rho1 = 3,
		P0 = 1,
		u0 = 1,
		v0 = 1,
	},
	--]]


	--initCond = 'sphere',

	--initCond = 'spiral',
	--initCondArgs = {torusGreaterRadius = .75, torusLesserRadius = .5},
	--[[ spiral with physically correct units, with 1 graph unit = 1 meter, and 1 simulation second = 1 second
	-- use this with euler-lingr + cylinder w/ r in [.5, 1], z in [-.25, .25]
	-- TODO instead of enforcing the boundary constraint, simulate across the full r=[0,1] domain and only impose a boundary on the fluid
	-- so that the GEM field can propagate outside the fluid
	initCondArgs = {
		solverVars={
			-- scale time down so that speedOfLight / units_m_per_s = 1 ... so second ~ 1/3e+8 ~ 3.33e-9
			second = 1 / constants.speedOfLight_in_m_per_s,

			speedOfLight = constants.speedOfLight_in_m_per_s, -- ... speedOfLight = 3e+8 / (1/(1/3e+8)) = 1
			divPsiWavespeed_g = constants.speedOfLight_in_m_per_s,
			divPhiWavespeed_g = constants.speedOfLight_in_m_per_s,

			-- 1 / units_m3_per_kg_s2 = second * second * kilogram / meter^3 ... so G ~ 6.67e-11 * 1e-17 ~ 6.67e-28
			gravitationalConstant = constants.gravitationalConstant_in_m3_per_kg_s2,

			-- change this so that rho=1 corresponds to the density of ... lead? mercury?
			kilogram = materials.Mercury.seaLevelDensity,
			kelvin = constants.K_to_C_offset + 12,	--materials.Mercury.boilingPoint,
		},

		-- with second=1, much higher than 1e-3*speedOfLight_in_m_per_s and we numerically explode
		-- so change second= higher to get it more stable
		v = constants.speedOfLight_in_m_per_s,
	},
	--]]

	--initCond = 'rarefaction wave',
	--initCond = 'Bessel',
	--initCond = 'jet',


	--initCond = 'Sod',
	--initCondArgs = {dim=cmdline.displayDim},
	--[[ real-world vars for Sod ... which are a few orders higher, and therefore screw up the backward-euler solver
	-- 		which means, todo, redo the backward euler error metric so it is independent of magnitude ... ?   seems I removed that for another numerical error reason.
	initCondArgs = {
		rhoL = 8 * materials.Air.seaLevelDensity,	-- kg / m^3
		PL = 10 * materials.Air.seaLevelPressure,	-- Pa = N / m^2 = kg / (m s^2)
		rhoR = materials.Air.seaLevelDensity,
		PR = materials.Air.seaLevelPressure,
		solverVars = {
			heatCapacityRatio = assert(materials.Air.heatCapacityRatio),
		},
	},
	--]]
	--[[
	some various initial conditions from 2012 Toro "The HLLC Riemann Solver" http://marian.fsik.cvut.cz/~bodnar/PragueSum_2012/Toro_2-HLLC-RiemannSolver.pdf
	Test	ρ L		u L		    p L		 ρ R		u R		    p R
	1		1.0		0.75		1.0		 0.125		0.0		    0.1
	2		1.0		-2.0		0.4		 1.0		2.0		    0.4
	3		1.0		0.0		    1000.0	 1.0		0.0		    0.01
	4		5.99924	19.5975		460.894	 5.99242	-6.19633	46.0950
	5		1.0		-19.59745	1000.0	 1.0		-19.59745	0.01
	6		1.4		0.0		    1.0		 1.0		0.0		    1.0
	7		1.4		0.1		    1.0		 1.0		0.1		    1.0
	--]]
	--[[ Sod vacuum test:
	initCondArgs = {
		rhoR = 0,
		PR = 0,
	},
	--]]


	initCond = 'rectangle',
	--initCond = 'Sedov',
	--initCond = 'Noh',
	--initCond = 'implosion',

	--initCond = 'Kelvin-Helmholtz',
	--[[
	initCondArgs = {
		noiseAmplitude = 1e-5,
		thickness = 1e-5,
	},
	--]]

	--initCond = 'Rayleigh-Taylor',	--FIXME ... get initial / static hydro potential working
	--initCond = 'Taylor-Green',	-- should only work with viscosity
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
	--initCond = 'self-gravitation - NGC 1560',	-- TODO still needs velocity
	--initCond = 'self-gravitation - NGC 3198',
	--initCond = 'self-gravitation test 1',
	--initCond = 'self-gravitation test 1 spinning',
	--initCond = 'self-gravitation test 2',		--FIXME
	--initCond = 'self-gravitation test 2 orbiting',
	--initCond = 'self-gravitation test 4',
	--initCond = 'self-gravitation soup',
	--initCond = 'self-gravitation Jeans, right?',


	--initCond = 'shallow water constant',
	--initCond = 'shallow water problem A',	-- boundary: v = reflect, h = freeflow
	--initCond = 'shallow water problem B',	-- boundary: v = reflect, h = freeflow
	--initCondArgs = {phi0 = 1},
	--initCond = 'shallow water problem C',	-- boundary = freeflow
	--TODO initCond = 'shallow water problem D',	-- lhs boundary = fixed, rhs boundary = freeflow
	--initCond = 'shallow water parabola',
	--initCond = '2003 Rogers',


	-- those designed for SRHD / GRHD from Marti & Muller 1998:
	--initCond = 'relativistic shock reflection',			-- FIXME.  these initial conditions are constant =P
	--initCond = 'relativistic blast wave test problem 1',
	--initCond = 'relativistic blast wave test problem 2',
	--initCond = 'relativistic blast wave interaction',		-- in 2D this only works with no limiter / lots of dissipation



	-- Maxwell: (doesn't work when run for euler, due to eqn env requirements)
	--initCond = 'Maxwell default',	-- arbitrary test to make sure it doesn't explode
	--initCond = 'Maxwell constant',
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

	-- initConds for twofluid-emhd stored in hydro/init/twofluid-emhd.lua:
	--initCond = 'two-fluid Brio-Wu', eqnArgs = {useEulerInitState=false},
	--initCond = 'GEM challenge', eqnArgs = {useEulerInitState=false},

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

	--initCond = 'black hole - isotropic - stuffed',

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
	--initCond = 'testbed - gauge wave',
	--initCond = 'testbed - gauge wave - diagonal',
	--initCond = 'testbed - linear wave',
	--initCond = 'testbed - linear wave - diagonal',
	--initCond = 'testbed - Gowdy',


	-- NLS
	--initCond = 'Gaussian',
	--initCond = 'Ring',
	--initCond = 'Oscillatory',
	---- piggybacking for my wave-finite-difference here: ----
	--initCond = 'Wave-FD Gaussian',
	--initCond = 'Wave-FD Bessel',

	-- multi-devices
	multiSlices = {cmdline.multiSlices or 2, 1, 1},
}
-- apply any cmdline 'initCondArgs' on top of the args.initCondArgs
if cmdline.initCondArgs or args.initCondArgs then
	args.initCondArgs = table(cmdline.initCondArgs, args.initCondArgs)
end

if cmdline.solver then self.solvers:insert(require('hydro.solver.'..cmdline.solver)(table(args, cmdline.solverArgs))) return end


-- simple wave equation, no time/space coupling via background metric


--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='wave'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='wave'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='rusanov', eqn='wave'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='wave', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='wave', wenoMethod='2008 Borges', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='wave', wenoMethod='2010 Shen Zha', order=5})))
--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='wave'})))


-- wave eqation with background spacetime metric
-- TODO explodes when addSource is used.


--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='wave_metric', eqnArgs={beta={'-y / (r * r)','x / (r * r)','0'}}, wenoMethod='1996 Jiang Shu', order=5})))

--[[ Acoustic black hole.  for use with cylinder grid?
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
self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='wave_metric'})))
--]]




-- shallow water equations


--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='shallow-water'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='shallow-water'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='shallow-water', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='shallow-water'})))

--[[ keep my configuration for the 2D world shallow water simulation in one place
-- roe seems to have less waves than hll
-- RK2 TVD seems to dissipate waves more than FE.  RK4 TVD even better.
-- and the bathymetry file is 4320x2160
self.solvers:insert(require 'hydro.solver.weno'(table(args, {
	eqn = 'shallow-water',
	dim = 2,
	useBathymetry = true,	-- TODO point to filename? right now filename is hardcoded into hydro/eqn/shallow-water.lua
	integrator = 'Runge-Kutta 4, TVD',
	wenoMethod = '1996 Jiang Shu',
	order = 5,
	coord = 'cartesian',
	coordArgs = {vectorComponent='cartesian'},
	mins = {-180, -90, -1},
	maxs = {180, 90, 1},
	gridSize = {432, 216},
})))
--]]


-- compressible Euler equations


--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler'})))

--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='euler', hllCalcWaveMethod=0})))	-- 'Davis direct'
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='euler', hllCalcWaveMethod=1})))	-- 'Davis direct bounded' -- this is the default hllCalcWaveMethod

--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='rusanov', eqn='euler'})))

--self.solvers:insert(require 'hydro.solver.fdsolver'(table(args, {eqn='euler'})))

--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='euler-hllc', eqn='euler', fluxArgs={hllcMethod=0}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='euler-hllc', eqn='euler', fluxArgs={hllcMethod=1}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='euler-hllc', eqn='euler', fluxArgs={hllcMethod=2}})))

-- NOTICE, these are very accurate with RK4, etc., but incur oscillations with Forward-Euler
-- TODO weno doesn't seem to work with self-gravitation
self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=5})))
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

-- testing different flux methods
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=5, fluxMethod='Lax-Friedrichs'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=5, fluxMethod='Roe'})))
-- FIXME: Marquina is broken:
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='1996 Jiang Shu', order=5, fluxMethod='Marquina'})))

--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', wenoMethod='2010 Shen Zha', order=5, fluxMethod='Roe'})))

-- for 11th WENO (2010 Shen Zha) once we reduce size below 6,6 it breaks
-- so TODO something about boundary conditions on WENO or something ... maybe an error
-- other than weno, this works fine with finite volume codes
--gridSize={64,1,1},


-- incompressible
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler', eqnArgs={incompressible=true}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='euler', eqnArgs={incompressible=true}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='euler-hllc', eqn='euler', eqnArgs={incompressible=true}})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='euler', eqnArgs={incompressible=true}, wenoMethod='2010 Shen Zha', order=5})))

-- viscosity in source terms as finite-different explicit update
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler', eqnArgs={viscosity='rhs-explicit'}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler', eqnArgs={viscosity='rhs-implicit'}})))
-- not yet finished:
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler', eqnArgs={viscosity='flux'}})))

-- incompressible + viscosity
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler', eqnArgs={incompressible=true, viscosity='rhs-explicit'}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler', eqnArgs={incompressible=true, viscosity='rhs-implicit'}})))

--[[ messing with viscosity settings
self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {
	flux = 'roe',
	eqn = 'euler',
	eqnArgs = {
		--incompressible = true,
		viscosity = 'rhs-explicit',
		shearViscosity = .001,
		heatConductivity = .02,
	}
})))
--]]

-- compressible Euler equations, based on primitive vars
-- only works with eigensystem-based solvers (since idk how to solve the flux vector of the primitive system)

--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler_prim'})))


-- Navier-Stokes-Wilcox:
-- TODO FIXME i think the eigenvector code isn't correct
-- my math & mathematica agree that the flux jacobian eigensystem is defective, therefore I don't know if a flux vector split method can be used ...
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='navstokes-wilcox'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='navstokes-wilcox', eqnArgs={incompressible=true}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='navstokes-wilcox'})))


-- compressible Euler equations - Burgers solver


-- f.e. and b.e. are working, but none of the r.k. integrators
-- PLM isn't implemented yet
-- neither is source term / poisson stuff
--self.solvers:insert(require 'hydro.solver.euler-burgers'(args))


-- compressible Euler fluid equations + de-Donder gauge linearized GR


--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler-lingr'})))


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

-- srhd incompressible.  (TODO explodes ... under which I.C? works fine with Sod.)
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='srhd', eqnArgs={incompressible=true}})))


-- general relativistic compressible hydrodynamics

-- TODO TODO this is a long way out of date
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
-- hmm, something is wrong, E waves propagating much faster than B waves, esp compared to glm-maxwell which looks good.

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
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='twofluid-emhd'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='twofluid-emhd', wenoMethod='1996 Jiang Shu', order=9})))	-- exploded...
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='twofluid-emhd', wenoMethod='2010 Shen Zha', order=5})))


-- here's another one: two-fluid emhd with de Donder gauge linearized general relativity
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='twofluid-emhd-lingr'})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='twofluid-emhd-lingr', wenoMethod='1996 Jiang Shu', order=5})))
--self.solvers:insert(require 'hydro.solver.weno'(table(args, {eqn='twofluid-emhd-lingr', wenoMethod='2010 Shen Zha', order=7})))


-- de Donder gauge + weak field limit + cartesian components, as a wave equation (+ TODO fluid as well)
-- TODO RENAME.  the twofluid lingr is gravitoelectromagnetism, which is a ,tt=0 approximation
--  while this is not.  it is a wave equation.
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='lingr'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='lingr'})))


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
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='z4_2008yano'})))

-- FIXME ? or not?  if I run 2D Alcubierre with this with cl-cpu then there's no problem.  Intel OpenCL, what's up?
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='z4'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='z4', eqnArgs={useShift='GammaDriverHyperbolic'}})))

--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='adm1d_v1'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='adm1d_v2'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='adm3d'})))	-- with 1D gaussian perturbation: explodes.  with 2D Alcubierre on cl-cpu 4 cores, looks like mem align errors
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='adm3d', eqnArgs={noZeroRowsInFlux=false}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='adm3d', eqnArgs={useShift='HarmonicShiftCondition-FiniteDifference'}})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='z4'})))
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='z4', eqnArgs={useShift='MinimalDistortionElliptic'}})))		-- with 1D gaussian perturbation: compile error
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='z4', eqnArgs={useShift='MinimalDistortionParabolic'}})))	-- with 1D gaussian perturbation: compile error
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='z4', eqnArgs={useShift='MinimalDistortionHyperbolic'}})))	-- with 1D gaussian perturbation: explodes
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='z4', eqnArgs={useShift='GammaDriverParabolic'}})))			-- with 1D gaussian perturbation: compile error
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='z4', eqnArgs={useShift='GammaDriverHyperbolic'}})))			-- with 1D gaussian perturbation: explodes

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

--self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn='bssnok-fd-senr'})))	-- default shift is HyperbolicGammaDriver

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
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='p2dfmt', meshfile='n0012_113-33.p2dfmt'}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2dcbrt', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2dcubed', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2dtwist', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='cylinder2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='cube3d', size={16, 16, 16}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='cylinder3d', size={16, 16, 16}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='sphere3d', size={16, 16, 16}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='torus3d', size={16, 16, 16}}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='image2d', extrude=1, image='blueprints/blueprint.png'}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', eqnArgs={incompressible=true}, mesh={type='image2d', extrude=1, image='blueprints/blueprint.png'}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='hll', eqn='euler', mesh={type='quad2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='euler-hllc', eqn='euler', mesh={type='quad2d', size={64, 64}}})))

-- TODO add boundary classes:

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='hll', eqn='euler', mesh={type='quad2d_with_cylinder_removed', size={32, 32}}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2d_with_cylinder_removed', size={32, 32}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2d_with_cylinder_removed', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2d_with_cylinder_removed', size={128, 128}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2d_with_cylinder_removed', size={256, 256}}})))

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', mesh={type='quad2d_with_box_removed', size={64, 64}}})))

-- TODO unstructured+incompressible
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler', eqnArgs={incompressible=true}, mesh={type='quad2d_with_cylinder_removed', size={64, 64}}})))

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

--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler-lingr', mesh={type='polar2d', size={64, 64}}})))
--self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {flux='roe', eqn='euler-lingr', mesh={type='cylinder3d', size={8, 8, 8}, mins={.5, 0, -.25}, maxs={1, 1, .25}}})))


if cmdline.naca then
	-- [=[ 1.25 degree angle of attack, mach 0.8, sea level pressure and density
	-- might be trying to reproduce the "I Do Like CFD" OssanWorld.com edu2d "case_steady_airfoil"
	--local theta = 0					-- cylinder uses 0
	local theta = math.rad(1.25)	-- naca airfoil uses 1.25
	--local machSpeed = 0.3			-- cylinder uses 0.3
	local machSpeed = 0.8			-- naca airfoil uses 0.8
	--local machSpeed = 0.95
	--local machSpeed = 1
	--local machSpeed = 1.2
	--local machSpeed = 2
	local kg = 1
	local m = 1
	local s = 1
	--local s = 10
	--local s = 100
	-- [[ using real units ...
	local Air = materials.Air
	local gamma = materials.Air.heatCapacityRatio
	local rho0 = Air.seaLevelDensity / (kg/m^3)		-- 1.2754 kg/m^3
	local v0 = Air.speedOfSound * machSpeed / (m/s)
	local P0 = Air.seaLevelPressure / (kg/(m*s^2))	-- 101325 (Pa = kg/(m s^2))
	--]]
	--[[ used in edu2d
	-- edu2d_module_ccfv_data_soln.f90, line 189 ... why is M_inf=0.8 used directly, and not M_inf * sqrt(gamma * P / rho) ?
	local gamma = 1.4
	local rho0 = 1
	local v0 = machSpeed / (m/s)
	local P0 = 1 / gamma
	--]]
	self.solvers:insert(require 'hydro.solver.meshsolver'(table(args, {
		flux = 'roe',
		--flux = 'hll',
		fluxArgs = {useEntropyFluxFix = true},
		fluxLimiter = 'donor cell',
		eqn = cmdline.eqn or 'euler',
		mesh = {
			--[[
			type = 'p2dfmt',	-- hmm, I don't know that p2dfmt has distinct boundary classes
			meshfile = 'n0012_113-33.p2dfmt',
			--]]
			-- [[
			type = 'edu2dgrid',
			--meshfile = 'cylinder.grid',
			meshfile = 'airfoil.grid',
			--]]
			boundary = {
				{	-- slip_wall
					name='mirror',
					args={restitution=0},
				},
				{	-- freestream
					name = 'fixed',
					args = {
						fixedCode = function(self, args, dst)
							local solver = args.solver
							return solver.eqn:template([[
{
	<?=prim_t?> W;	//force entire state to initial state
	W.rho = <?=rho0?>;
	W.P = <?=P0?>;
	W.v.x = <?=math.cos(theta) * v0?>;
	W.v.y = <?=math.sin(theta) * v0?>;
	W.v.z = 0;
	W.ePot = 0;
	<?=consFromPrim?>(<?=dst?>, solver, &W, <?=face?>->pos);
}
]],
								table(args, {
									rho0 = rho0,
									v0 = v0,
									P0 = P0,
									theta = theta,
								}):setmetatable(nil)
							),
							-- 2nd return is list of required modules
							{assert(solver.eqn.symbols.consFromPrim)}
						end,
					}
				},
			},
			{
				-- "outflow_subsonic", used with the cylinder test
				name = 'fixed',
				args = {
					fixedCode = function(self, args, dst)
						local solver = args.solver
						return solver.eqn:template([[
{
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, <?=src?>, <?=face?>->pos);
	W.P = <?=P0?>;	//force pressure to initial pressure
	<?=consFromPrim?>(<?=dst?>, solver, &W, <?=face?>->pos);
}
]],
							table(args, {P0=P0}):setmetatable(nil)
						),
						-- 2nd return is list of required modules
						{
							assert(solver.eqn.symbols.consFromPrim),
							solver.eqn.symbols.primFromCons,
						}
					end,
				},
			},
		},
		initCond = 'constant',
		initCondArgs = {
			solverVars = {
				meter = m,
				second = s,
				kilogram = kg,
				heatCapacityRatio = gamma,
			},
			rho = rho0,
			v = {
				math.cos(theta) * v0,
				math.sin(theta) * v0,
				0,
			},
			P = P0,
			-- so if we want to reduce this to 1, we can scale down seconds ...
		},
		-- [[
		--cfl = 1,
		--cfl = 0.9,	-- I think this is what "I do like CFD" uses ... maybe ... they have an extra 0.5 in the denom tho, makes me think 1.8 is their cfl ... and they also evaluate it for "steady" on a per-cell basis
		--cfl = 1.8,	-- the "/ 0.5" in "I do like cfd" ... wouldn't that be the same as x2? not exactly...
		--cfl = .45,
		-- the first unsteady dt of airfoil in edu2d is 2.8647764373497124E-004 (with their CFL=0.9) ... let's see what cfl I need to use with my implementation of face-based cfl to reproduce this ...
		-- with my implementation (without a 0.5 in the denom ... is that there because dim = 2?  i thought it was divide-by-dim, not multiply-by-dim ...) and cfl=1 we get a dt=0.00036828253563196
		-- that means to get the same dt, we should use cfl=0.77787463704567
		cfl = 0.77787463704567,	-- produces dt=0.00028647764373497
		--integrator = 'forward Euler',
		integrator = 'Runge-Kutta 2, TVD',
		--integrator = 'Runge-Kutta 4',
		--]]
		--[[	-- not yet working with meshsolver
		cfl = 4,
		--integrator = 'backward Euler',
		--integrator = 'backward Euler, CPU',
		--]]
	})))
	--]=]
end

-- NEXT BIG TODO
-- * make meshsolver and gridsolver separate options
-- * make grid coordinate chart separate of vector component coordinate chart
-- so in the end, the user can choose the geometry: mesh (w/coordinate mapping) vs grid
--		The only dif between grid and mesh solver is that grid has easy n'th order finite difference stencils, thanks to its parameterization
-- 		I can unify this with mesh if I store per-cell the d/dx^i(x + j*dx^i) in the i'th face direction to the j'th step.
--		Then I could rewrite the current higher order PLM/WENO stuff for the MeshSolver
-- in fact, the biggest difference between mesh and grid solver is the n'th spatial derivative calculations, so if this can be abstracted then the two can be combined


-- multi GPU
-- how about 'composite grid' instead of 'chopped up'?
--self.solvers:insert(require 'hydro.solver.choppedup'(table(args, {flux='roe', eqn='euler', subsolverClass=require 'hydro.solver.fvsolver'})))


-- composite equations.  better than composite solver. less kernel calls.
-- single with distinct cons_t & prim_t
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='composite', eqnArgs={subeqns={'euler'}}})))
-- single with matching cons_t & prim_t
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='composite', eqnArgs={subeqns={'maxwell'}}})))
-- multiple w/ distinct cons_t & prim_t
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='composite', eqnArgs={subeqns={'euler', 'euler'}}})))
-- multiple w/ matching cons_t & prim_t
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='composite', eqnArgs={subeqns={'maxwell', 'maxwell'}}})))
-- multiple w/ separate, distinct & matching cons_t & prim_t
-- TODO FIXME getting memory alignment errors, even though the structs in 'checkStructSizes' say they line up.
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='composite', eqnArgs={subeqns={'euler', 'maxwell'}}})))
-- two fluid plasma eventually
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='composite', eqnArgs={subeqns={'euler', 'euler', 'maxwell'}}})))


-- the start of AMR
--self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler'})))
--self.solvers:insert(require 'hydro.solver.amr'(require 'hydro.solver.fvsolver')(table(args, {flux='roe', eqn='euler'})))


-- [=[ reproducing 2009 Alic, Bona, Bona-Casas"Towards a gauge-polyvalent numerical relativity code"
local dim = cmdline.dim or 1
local args = {
	app = self,

	--eqn = 'adm3d',
	eqn = 'z4',

	integrator = 'forward Euler',
	--integrator = 'Runge-Kutta 3, TVD',	-- p.20, eqn B.1
	--integrator = 'Runge-Kutta 4',
	dim = dim,
	cfl = .5/dim,			-- no mention of cfl or timestep ...
	fluxLimiter = cmdline.fluxLimiter or 'donor cell',	-- I didn't read much into what kind of flux/slope limiter was used
	--fluxLimiter = cmdline.fluxLimiter or 'minmod',
	--fluxLimiter = cmdline.fluxLimiter or 'Lax-Wendroff',
	--fluxLimiter = cmdline.fluxLimiter or 'superbee',

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
	mins = {.1, 0, 0},
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
	--[[ sphere sinh radial with parameters from 2009 Alic (right?)
	coord = 'sphere_sinh_radial',
	coordArgs = {
		-- TODO sort this out
		-- TODO do I have mem write / unwritten vars in "holonomic"?  cuz there seem to be errors that persist past reset()
		-- TODO move cell_area and cell_volume calcs into cell_t fields
		--vectorComponent = 'cartesian',
		vectorComponent = 'anholonomic',	-- ... these settings also influence the finite volume area/volume calculations (in terms of the vector components) ...
		--vectorComponent = 'holonomic',	-- our tensor components are holonomic ... except the partial / 1st order state variables, like a_k, d_kij

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
	-- [[ sphere_sinh_radial but with SENR parameters ... for SENR init conds
	coord = 'sphere_sinh_radial',
	coordArgs = {
		--vectorComponent = 'cartesian',
		vectorComponent = 'anholonomic',
		--vectorComponent = 'holonomic',
		--sinh_w = .15,
		sinh_w = math.sinh(0.0916845),	-- 2017 Ruchlin, Fig 4, w = 0.0916845
		amplitude = 1000,				-- 2017 Ruchlin, Fig 4, rmax / M = 1000
	},
	mins = {0, 0, 0},
	maxs = {1, math.pi, 2*math.pi,},
	gridSize = cmdline.gridSize or ({
		{200, 1, 1},	-- 2017 Ruchlin, Fig 4, Nx1 = 200, Nx2 = 2, Nx3 = 2 ... but I'm not bothering with the finite-difference min width of 2 stuff
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

	--initCond = 'Minkowski',				-- stable with 2009Alic-z4 (with derivative vars initd to zero)
	--initCond = 'SENR Minkowski',			-- stable with 2009Alic-z4 (with derivative vars initd to zero)
	--initCond = 'plane gauge wave',
	initCond = 'SENR UIUC',				-- 2009Alic-z4 (with derivative vars initd to zero) coord=sphere_sinh_radial with SENR UIUC runs indefinitely.  coord=sphere runs until t=200 then e.h.? hits rhs and explodes.
	--initCond = 'UIUC',					-- but why does this one run so slow? smaller alphas means lower cfls? too close to zero?
	--initCond = 'SENR BrillLindquist',
	--initCond = 'black hole - Schwarzschild',
	--initCond = 'black hole - isotropic - stuffed',	-- TODO FIXME
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
	--[[
	initCond = 'gaussian perturbation',
	-- override these from above
	mins = {-1,-1,-1},
	maxs = {1,1,1},
	--]]
	--[[
	initCond = 'Alcubierre warp bubble',
	--initCondArgs = {R=.5, sigma=8, speed=.1},	-- sub-luminal
	--initCondArgs = {R=.5, sigma=8, speed=1.1},		-- super-luminal 1.1x
	-- override these from above
	mins = {-1,-1,-1},
	maxs = {1,1,1},
	-- quadratic is having numerical problems at the borders
	boundary = {
		xmin = 'freeflow',
		xmax = 'freeflow',
		ymin = 'freeflow',
		ymax = 'freeflow',
		zmin = 'freeflow',
		zmax = 'freeflow',
	},
	gridSize = {64,64,1},
	--]]

	flux = 'hll',
	--flux = 'rusanov',
}
-- comparing hll solvers
if cmdline['2009Alic-adm'] then
	self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {eqn='adm3d'}, cmdline.solverArgs)))
end
if cmdline['2009Alic-z4'] then
	self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {eqn='z4', eqnArgs={useShift=cmdline.useShift}, initCond=cmdline.initCond}, cmdline.solverArgs)))
end
if cmdline['2009Alic-z4_2008yano'] then
	self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {eqn='z4_2008yano'}, cmdline.solverArgs)))
end
if cmdline['2009Alic-bssnok-fd-senr'] then
	self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn='bssnok-fd-senr'}, cmdline.solverArgs)))
end
if cmdline['2009Alic-bssnok-fd-num'] then
	self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn='bssnok-fd-num'}, cmdline.solverArgs)))
end
if cmdline['2009Alic-bssnok-fd-sym'] then
	self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn='bssnok-fd-sym'}, cmdline.solverArgs)))
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

		--cflMethod = cmdline.cflMethod or '2008 Alcubierre',
		--cflMethod = cmdline.cflMethod or '2013 Baumgarte et al, eqn 32',
		cflMethod = cmdline.cflMethod or '2017 Ruchlin et al, eqn 53',
	},
	dim = dim,

	integrator = cmdline.integrator or 'Runge-Kutta 4',
	--integrator = cmdline.integrator or 'backward Euler',
	--integrator = cmdline.integrator or 'backward Euler, CPU',	-- debugging.   seems that, for grid sizes too small, B.E. GPU fails.  i think because the reduce() gpu function isn't set up for lower bounds of buffer sizes.
	--integratorArgs = {verbose=true},
	cfl = cmdline.cfl or .5,

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
	coord = 'sphere_sinh_radial',
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
		--{32, 2, 2},		-- SENR sphere_sinh_radial uses this by default
		{80, 80, 2},		-- this works well for BrillLindquist sphere_sinh_radial when viewing the xz slice
		--{128, 2, 2},
		--{128, 32, 2},
		--{400, 64, 2},

		--{200, 2, 2},

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


	--initCond = 'Minkowski',	-- TODO sphere_sinh_radial

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
	--initCond = 'SENR UIUC',					-- single black hole. bssnok-fd-num explodes because H diverges at t=13 ... when partial_phi_l diverges at the same rate ... because of its r=0 value?
	initCond = 'SENR BrillLindquist',			-- two merging head-on.
	--initCond = 'SENR BoostedSchwarzschild',
	--initCond = 'SENR StaticTrumpet',

	-- multi-devices
	multiSlices = {cmdline.multiSlices or 3, 1, 1},
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
	self.solvers:insert(require 'hydro.solver.bssnok-fd'(table(args, {eqn = 'bssnok-fd-senr', initCond = cmdline.initCond})))
end
if cmdline['bssnok-fd-senr-multi'] then
	self.solvers:insert(require 'hydro.solver.choppedup'(table(args, {eqn = 'bssnok-fd-senr', subsolverClass = require 'hydro.solver.bssnok-fd'})))
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
if cmdline['z4-hll'] then
	self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='hll', eqn='z4', eqnArgs = {useShift = cmdline.useShift}, initCond = cmdline.initCond})))
end
--]=]

-- [=[ ok here's my 2D z4 finite volume brill lindquist solver
if cmdline.z4bl then
	local dim = 2
	local args = {
		app = self,
		integrator = 'forward Euler',
		dim = dim,
		cfl = .5/dim,
		fluxLimiter = cmdline.fluxLimiter or 'donor cell',
		--fluxLimiter = cmdline.fluxLimiter or 'minmod',
		--fluxLimiter = cmdline.fluxLimiter or 'Lax-Wendroff',
		--fluxLimiter = cmdline.fluxLimiter or 'superbee',

		coord = 'sphere_sinh_radial',
		coordArgs = {
			vectorComponent = 'anholonomic',	-- ... these settings also influence the finite volume area/volume calculations (in terms of the vector components) ...
			amplitude = 1000,
			sinh_w = .15
		},
		mins = {0, 0, 0},
		maxs = {
			1,
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

		initCond = 'SENR BrillLindquist',
		flux = 'hll',
		eqn = 'z4',
		eqnArgs = {useShift = 'GammaDriverHyperbolic'},
	}
	self.solvers:insert(require 'hydro.solver.fvsolver'(args))
end
--]=]

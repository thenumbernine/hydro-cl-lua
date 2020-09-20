-- slowly making this the home of all initial conditions...
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local math = require 'ext.math'
local clnumber = require 'cl.obj.number'
local template = require 'template'
local materials = require 'hydro.materials'
local InitCond = require 'hydro.init.init'
local real = require 'hydro.real'
local constants = require 'hydro.constants'

local common = require 'hydro.common'
local xNames = common.xNames
local minmaxs = common.minmaxs

local function RiemannProblem(initCond)
	local WL, WR = initCond[1], initCond[2]

	local function getInitCondFieldName(i,varname)
		local suffix = ({'L','R'})[i]
		return varname:gsub('%.','')..suffix
	end

	initCond.guiVars = initCond.guiVars or {}
	for i,WLR in ipairs(initCond) do
		for name, value in pairs(WLR) do
			table.insert(initCond.guiVars, {
				name = getInitCondFieldName(i,name),
				value = value,
			})
		end
	end
	
	initCond.init = function(self, solver, args)
		if args then
			self.overrideDim = args.dim
		end
	end
	
	initCond.getInitCondCode = function(self, solver)
		local function build(i)
			return table.map(initCond[i], function(_,name,t)
				return '\t\t'..name..' = initCond->'..getInitCondFieldName(i,name)..';', #t+1
			end):concat'\n'
		end
		return template([[
	
	bool lhsSod = true<?
for i=1,overrideDim or solver.dim do
	local xi = xNames[i]
?> && x.<?=xi?> < mids.<?=xi?><?
end
?>;

	if (lhsSod) {
<?=build(1)?>
	} else {
<?=build(2)?>
	}
]], 		{
				solver = solver,
				build = build,
				xNames = xNames,
				overrideDim = self.overrideDim,
			})
	end


	-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
	-- TODO it might be more efficient to hand this function the array, and have it return an array
	-- or at least have it calculate certain values up front before iterating across all x's
	initCond.exactSolution = function(solver, x, t)
		local solverPtr = solver.solverPtr
		local initCondPtr = solver.initCondPtr
		-- TODO initial condition object to share these values with initialization
		local rhoL = initCondPtr.rhoL
		local rhoR = initCondPtr.rhoR
		local PL = initCondPtr.PL
		local PR = initCondPtr.PR
		local vL = 0
		local vR = 0
		local gamma = solverPtr.heatCapacityRatio
		
		local muSq = (gamma - 1)/(gamma + 1)
		local K = PL / rhoL^gamma

		local CsL = math.sqrt(gamma * PL / rhoL)
		local CsR = math.sqrt(gamma * PR / rhoR)

		local solveP3 = function()
			-- hmm, for some reason, using the symmath-compiled code is resulting in memory corruption
			-- it's not like symmath compile performs any ffi allocations ... maybe something is getting freed prematurely?
			-- or maybe I'm passing in a cdata that is a number, and luajit is not agreeing with some implicit conversion to a Lua number somewhere?  I don't know...
			-- so to fix this, I'll just print out the code and inline it myself
			--[=[ using symmath-compiled functions
			local symmath = require 'symmath'
			local P3, PL, PR, CsL, CsR, gamma = symmath.vars('P3', 'PL', 'PR', 'CsL', 'CsR', 'gamma')
			local f = -2*CsL*(1 - (P3/PL)^((-1 + gamma)/(2*gamma)))/(CsR*(-1 + gamma)) + (-1 + P3/PR)*((1 - muSq)/(gamma*(muSq + P3/PR)))^.5
			local df_dP3 = f:diff(P3)()	
			local vars = {P3, PL, PR, CsL, CsR, gamma}
			local f_func, f_code = f:compile(vars)
			local df_dP3_func, df_dP3_code = df_dP3:compile(vars)
			--print(f_code)
			--print(df_dP3_code)
			--os.exit()
			--]=]
			local P3 = .5 * (PL + PR)
			local epsilon = 1e-16	-- this is the limit for the sod.f before it oscillates
			while true do
				--[=[ using symmath-compiled functions
				local dP3 = -f_func(P3) / df_dP3_func(P3)
				--]=]
				-- [=[ using inlining
				local f = ((((-2 * CsL) * (1 - ((P3 / PL) ^ ((-1 + gamma) / (2 * gamma))))) / (CsR * (-1 + gamma))) + ((-1 + (P3 / PR)) * ((0.75 / (gamma * (0.25 + (P3 / PR)))) ^ 0.5))) 
				local df_dP3 = ((-((((((1.5 * math.sqrt(0.75) * CsR * PR * (gamma ^ 1.5)) - ((0.75 ^ 1.5) * CsR * PR * math.sqrt(gamma))) - ((0.75 ^ 1.5) * CsR * (gamma ^ 2.5) * PR)) - (0.5 * P3 * math.sqrt(0.75) * CsR * math.sqrt(gamma))) - (0.5 * P3 * math.sqrt(0.75) * CsR * (gamma ^ 2.5))) + (((P3 * math.sqrt(0.75) * CsR * (gamma ^ 1.5)) - (0.25 * (PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ (((-1) - gamma) / (2 * gamma))) * (PR ^ 1.5) * math.sqrt((P3 + (0.25 * PR))))) - ((PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ ((-(1 - gamma)) / (2 * gamma))) * math.sqrt(PR) * math.sqrt((P3 + (0.25 * PR))))) + (0.5 * (PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ (((-1) - gamma) / (2 * gamma))) * (PR ^ 1.5) * gamma * math.sqrt((P3 + (0.25 * PR)))) + (((2 * (PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ ((-(1 - gamma)) / (2 * gamma))) * math.sqrt(PR) * gamma * math.sqrt((P3 + (0.25 * PR)))) - (0.25 * (PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ (((-1) - gamma) / (2 * gamma))) * (gamma ^ 2) * (PR ^ 1.5) * math.sqrt((P3 + (0.25 * PR))))) - ((PL ^ ((1 - gamma) / (2 * gamma))) * CsL * (P3 ^ ((-(1 - gamma)) / (2 * gamma))) * (gamma ^ 2) * math.sqrt(PR) * math.sqrt((P3 + (0.25 * PR))))))) / (math.sqrt(PR) * CsR * ((P3 + (0.25 * PR)) ^ 1.5) * gamma * ((1 - (2 * gamma)) + (gamma ^ 2)))) 
				local dP3 = -f / df_dP3
				--]=]
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

	
		-- between regions 1 and 2
		local s1 = -CsL	

		-- between regions 2 and 3
		-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
		local s2 = -vtail

		local s3 = v3	-- between regions 3 and 4

		-- between regions 4 and 5 ...
		local s4 = vshock

		--print('wavespeeds:',s1,s2,s3,s4)

		local rho, vx, P
		local xi = x / t
		if xi < s1 then
			rho = rhoL
			vx = vL
			P = PL
		elseif xi < s2 then
			vx = (1 - muSq) * (x/t + CsL)
			
			-- Dullemon:
			--rho = (rhoL^gamma / (gamma * PL) * (v2(x) - x/t)^2)^(1/(gamma-1))
			-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
			rho = rhoL * (-muSq * (x / (CsL * t)) + (1 - muSq))^(2/(gamma-1))

			-- Dullemon:
			--P = K * rho2^gamma
			-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
			P = PL * (-muSq * (x / (CsL * t)) + (1 - muSq)) ^ (2*gamma/(gamma-1))
		elseif xi < s3 then
			rho = rho3
			vx = v3
			P = P3
		elseif xi < s4 then
			rho = rho4
			vx = v4
			P = P4
		else
			rho = rhoR
			vx = vR
			P = PR
		end

		local vy = 0
		local vz = 0
		local EInt = P / (gamma - 1)
		local EKin = .5 * rho * (vx*vx + vy*vy + vz*vz)
		local ETotal = EKin + EInt
		return rho, rho * vx, rho * vy, rho * vz, ETotal
	end

	return initCond
end

local function quadrantProblem(initCond)
	initCond.getInitCondCode = function(self, solver)
		-- [[ specific to 2002 Kurganov, Tadmor, "Solution of Two-Dimensional Riemann Problems for Gas Dynamics without Riemann Problem Solvers"
		solver.cfl = .475
		solver:setBoundaryMethods'freeflow'
		--]]
		local function build(i)
			local q = initCond[i]
			return table.map(q, function(v,k,t)
				return k..'='..v..';', #t+1
			end):concat' '
		end
		return template([[
	bool xp = x.x > mids.x;
	bool yp = x.y > mids.y;
	if (yp) {
		if (xp) {
			<?=build(1)?>
		} else {
			<?=build(2)?>
		}
	} else {
		if (!xp) {
			<?=build(3)?>
		} else {
			<?=build(4)?>
		}
	}
]],			{
				build = build,
			})
	end
	return initCond
end

-- right now 'center' is provided in cartesian coordinates (i.e. post-applying coordMap)
local SelfGravProblem = class()

function SelfGravProblem:getRadiusCode(source)
	return clnumber(source.radius)
end

function SelfGravProblem:init(args)
	self.args = args
	self.getRadiusCode = args.getRadiusCode
end

SelfGravProblem.depends = {'coordMap'}
function SelfGravProblem:getInitCondCode(initCond, solver)
	local args = self.args

	solver.useGravity = true

	return template([[
	rho = .1;
	P = 1;
	//notice about initializing random velocity -- it isn't uniform about .5 so it will pull left
	//v.x = .2 * (U->m.x - .5);	//U is initialized to random()
	//v.y = .2 * (U->m.y - .5);
	//v.z = .2 * (U->m.z - .5);
	
	<? for i,source in ipairs(sources) do ?>{
		real3 xc = coordMap(x);
		real3 delta = real3_sub(xc, _real3(
			<?=clnumber(source.center[1])?>,
			<?=clnumber(source.center[2])?>,
			<?=clnumber(source.center[3])?>));
		real distSq = real3_lenSq(delta);
		real radius = <?=self:getRadiusCode(source)?>;
		if (distSq < radius * radius) {
			<?=source.inside or 'rho = P = 1;'?>
		}
	}<? end ?>
]], {
		self = self,
		sources = args.sources,
		clnumber = clnumber,
	})
end

--[[
args:
	solver = solver
	side = which side, default xmin
	dir = which direction the wave will point. default 'x'
	amplitude = wave amplitude.  default 1
	frequency = wave frequency.  default 1
--]]
local function addMaxwellOscillatingBoundary(args)
	local solver = assert(args.solver)
	local side = args.side or 'xmin'
	local amplitude = args.amplitude or 1
	local frequency = args.frequency or 1
	local dir = args.dir or 'y'

	-- TODO addBoundaryOption?  only do this once?
	local BoundaryOscillating = class(solver.Boundary)
	BoundaryOscillating.name = 'oscillating'
	function BoundaryOscillating:getCode(args)
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		local U = 'buf['..args.index(
			side:sub(-3) == 'min' and 'j+1' or (gridSizeSide..'-numGhost-1-j')
		)..']'
		
		-- TODO put the old code here
		if not args.fields then
			return 
				--oldxmin(args) .. 
				template([[
<?
	local epsSrc = 
		(require 'hydro.eqn.twofluid-emhd-lingr'.is(eqn) 
			or require 'hydro.eqn.glm-maxwell'.is(eqn)
			or require 'hydro.eqn.maxwell'.is(eqn)
		)
		and U or 'solver'
?>
	
	<?=U?>.B = <?=vec3?>_zero;
	<?=U?>.D = <?=vec3?>_zero;
	<?=U?>.D.<?=dir?> = <?=real_mul?>(
		<?=mul?>(<?=epsSrc?>.sqrt_1_eps, <?=epsSrc?>.sqrt_1_eps),
		<?=amplitude?> * sin(2. * M_PI * <?=frequency?> * t)); 
]], 		table(solver.eqn:getEnv(), {
				U = U,
				eqn = self.eqn,
				frequency = clnumber(frequency),
				amplitude = clnumber(amplitude),
				dir = dir,
			}))
		else
			return ''
		end
	end
	
	-- this args is only for the UBuf boundary program -- not calle for the Poisson boundary program
	local oldGetBoundaryProgramArgs = solver.getBoundaryProgramArgs
	function solver:getBoundaryProgramArgs()
		-- i'm completely overriding this
		-- so I can completely override boundaryMethods for the solver boundary kernel
		-- yet not for the poisson boundary kernel
--		local boundaryMethods = table(self.boundaryMethods)
		-- TODO get the oscillations on 2D 256x256 in the upper left corner to stop
		--local oldxmin = select(2, next(solver.boundaryOptions[boundaryMethods.xmin]))
		self.boundaryMethods[side] = BoundaryOscillating()
		
		local args = oldGetBoundaryProgramArgs(self)
		
		-- same as super 
		-- except with extraAgs
		-- and using boundaryMethods instead of self.boundaryMethods
		-- (which I leave alone so Poisson can deal with it)
		
		-- TODO just give the solver a parameter 't' ?
		-- give it a parameter 'dt' while you're at it
		-- that would save on a few kernel parameters
		args = table(args)
		args.extraArgs = {'real t'}

		return args
	end
	
	-- this runs before refreshBoundaryProgram, so lets hijack refreshBoundaryProgram and add in our time-based boundary conditions
	local oldBoundary = solver.boundary
	function solver:boundary()
assert(self.t)
		for _,obj in ipairs(self.boundaryKernelObjs) do
			obj.obj:setArg(2, real(self.t))
		end
		oldBoundary(self)
	end
end

local initConds = table{
	{
		name = 'constant',
		guiVars = {
			{name='rho0', value=1},
			{name='vx0', value=0},
			{name='vy0', value=0},
			{name='vz0', value=0},
			{name='P0', value=1},
		},	
		solverVars = {
			heatCapacityRatio = 7/5,
		},
		init = function(self, solver, args)
			self.initCondArgs = args
		end,
		getInitCondCode = function(self, solver)
			local args = self.initCondArgs
			if args then
				local found
				if args.rho then self.guiVars.rho0.value = args.rho found = true end
				if args.P then self.guiVars.P0.value = args.P found = true end
				if args.v then 
					if args.v[1] then self.guiVars.vx0.value = args.v[1] found = true end
					if args.v[2] then self.guiVars.vy0.value = args.v[2] found = true end
					if args.v[3] then self.guiVars.vz0.value = args.v[3] found = true end
				end
			end	
		
			return [[
	rho = initCond->rho0;
	v.x = initCond->vx0;
	v.y = initCond->vy0;
	v.z = initCond->vz0;
	P = initCond->P0;
]]
		end,
	},
	{
		name = 'random',
		getInitCondCode = function(self, solver)
			solver.useGravity = true
			return [[
	rho = U->rho + 1.;
#if dim == 2
	v = _real3(
		U->m.x * cos(U->m.y * 2. * M_PI),
		U->m.x * sin(U->m.y * 2. * M_PI),
		2. * U->m.z - 1.
	);
#else
	// TODO spherical random distribution using asin and stuff
	v = real3_sub(U->m, _real3(.5, .5, .5));
#endif
	P = U->ETotal + 1.;
]]
		end,
	},
	{
		name = 'linear',
		solverVars = {
			heatCapacityRatio = 7/5,
		},
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			return [[
	real3 xc = coordMap(x);
	rho = 1. + xc.x;
	P = 1. + xc.x;
]]
		end,
	},
	
	-- 2017 Zingale section 7.9.3
	{
		name = 'gaussian',
		guiVars = {
			{name = 'rho0', value = 1e-3},
			{name = 'rho1', value = 1},
			{name = 'sigma', value = .1},
			{name = 'u0', value = 0},
			{name = 'v0', value = 0},
			{name = 'P0', value = 1e-6},
			{name = 'x0', value = -.5},
			{name = 'y0', value = -.5},
			{name = 'z0', value = 0},
		},
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			return template([[
	real3 xc = coordMap(x);
	//real xSq = real3_lenSq(xc);
	real xSq = real3_lenSq(real3_sub(xc, 
		_real3(
			initCond->x0,
			initCond->y0, 
			initCond->z0
		)));
	rho = (initCond->rho1 - initCond->rho0) * exp(-xSq / (initCond->sigma*initCond->sigma)) + initCond->rho0;
	v.x = initCond->u0;
	v.y = initCond->v0;
	P = initCond->P0;
]],		{
			clnumber = clnumber,
		})
		end,
	},

	{
		-- boundary waves seem to mess with this, 
		-- otherwise it looks like a wave equation solution
		name = 'Bessel',
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			return [[
	real r = coordMapR(x);
	real3 xc = coordMap(x);
	real3 n = real3_real_mul(xc, 1. / r);
	
	//for wave equation I'm mapping rho to Phi_,t = Pi
	rho = 0;

	//for wave equation I'm mapping v_i to Phi_,i = Psi_i
	//and I'm doing it in Cartesian coordinates for some reason
	// but I want the Bessel function to vary in radial coordinates
	// so ... I have to initialize this to d/dr J0(r)
	// and then rotate it to point along e_r
	// so, using -J1(r) = d/dr J0(r) according to Wiki (which says surpriringly little about Bessel function derivatives)
	real vr = -BESSJ1(r);
	v = real3_real_mul(n, vr);
]]
		end,
	},



	{
		name = 'advect wave',
		mins = {0,0,0},
		maxs = {1,1,1},
		solverVars = {
			heatCapacityRatio = 7/5,
		},
		guiVars = {
			{name = 'v0x', value = 1},
			{name = 'v0y', value = 0},
			{name = 'rho0', value = 1},
			{name = 'rho1', value = 3.2e-1},
			{name = 'P0', value = 1},
		},
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			solver:setBoundaryMethods{
				xmin = 'periodic',
				xmax = 'periodic',
				ymin = 'periodic',
				ymax = 'periodic',
				zmin = 'periodic',
				zmax = 'periodic',
			}
			return [[
	real3 xc = coordMap(x);
	real xmin = solver->mins.x;
	real xmax = solver->maxs.x;
	real width = xmax - xmin;
	real k0 = 2. * M_PI / width;
	rho = initCond->rho0 + initCond->rho1 * sin(k0 * (xc.x - xmin));
	v.x = initCond->v0x;
	v.y = initCond->v0y;
	v.z = 0;
	P = initCond->P0;
	ePot = 0;
]]
		end,
		-- TODO combine this with above, use a parser / transpiler to convert between Lua and OpenCL, and just write one equation?
		-- TODO TODO do this with all math everywhere, and analyze the dependency graph of variables and automatically slice out what GPU calculations should be buffered / automatically inline equations
		exactSolution = function(solver, x, t)
			local solverPtr = solver.solverPtr
			local initCondPtr = solver.initCondPtr
			local k0 = 2 * math.pi / (solverPtr.maxs.x - solverPtr.mins.x)
			local rho = initCondPtr.rho0 + initCondPtr.rho1 * math.sin(k0 * (x - t))
			local mx = 1 * rho
			local my = 0
			local mz = 0
			local P = 1
			-- hmm, only good with cartesian geometry
			local mSq = mx * mx + my * my + mz * mz 
			local EKin = .5 * mSq  / rho
			local EInt = P / (solverPtr.heatCapacityRatio - 1)
			local ETotal = EKin + EInt
			return rho, mx, my, mz, ETotal
		end,
	},

	-- test case vars
	RiemannProblem{
		name = 'Sod',
		-- L = high pressure / density
		{rho = 1, P = 1},
		-- R = low pressure / density
		{rho = .125, P = .1},
		solverVars = {
			heatCapacityRatio = 5/3,
		},
	},

-- [[ real-world vars ... which are a few orders higher, and therefore screw up the backward-euler solver
-- 		which means, todo, redo the backward euler error metric so it is independent of magnitude ... ?   seems I removed that for another numerical error reason.
	RiemannProblem{
		name = 'Sod with physical units',
		{
			rho = 8 * materials.Air.seaLevelDensity,	-- kg / m^3
			P = 10 * materials.Air.seaLevelPressure,	-- Pa = N / m^2 = kg / (m s^2)
		},
		{
			rho = materials.Air.seaLevelDensity,
			P = materials.Air.seaLevelPressure,
		},
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


	{	-- just like Brio-Wu, but centered instead of to one side
		name = 'rectangle',
		guiVars = {
			{name = 'rhoL', value = 1},
			{name = 'rhoR', value = .125},
			{name = 'PL', value = 1},
			{name = 'PR', value = .1},
			{name = 'BxL', value = .75},
			{name = 'BxR', value = .75},
			{name = 'ByL', value = 1},
			{name = 'ByR', value = -1},
			{name = 'BzL', value = 0},
			{name = 'BzR', value = 0},
		},
		solverVars = {
			heatCapacityRatio = 2,
		},
		getInitCondCode = function(self, solver)
			return template([[
	lhs = true 
<?
for i=1,solver.dim do
	local xi = xNames[i]
?> 		&& x.<?=xi?> > .75 * solver->mins.<?=xi?> + .25 * solver->maxs.<?=xi?>
		&& x.<?=xi?> < .25 * solver->mins.<?=xi?> + .75 * solver->maxs.<?=xi?>
<?
end
?>;
	
	rho = lhs ? initCond->rhoL : initCond->rhoR;
	P = lhs ? initCond->PL : initCond->PR;
	B.x = lhs ? initCond->BxL : initCond->BxR;
	B.y = lhs ? initCond->ByL : initCond->ByR;
	B.z = lhs ? initCond->BzL : initCond->BzR;
]], 	{
			solver = solver,
			xNames = xNames,
		}) 
		end,
	},

	{
		name = 'Sedov',
		guiVars = {
			{name = 'rho0', value = 1},
			{name = 'P0', value = 1},
			{name = 'P1', value = 1e+3},
		},
	
		getInitCondCode = function(self, solver)
			return [[
	rho = initCond->rho0;
	P = (i.x == solver->gridSize.x/2 && i.y == solver->gridSize.y/2 && i.z == solver->gridSize.z/2) ? initCond->P1 : initCond->P0;
]]
		end,
	},
	
	-- http://www-troja.fjfi.cvut.cz/~liska/CompareEuler/compare8/node36_ct.html
	{
		name = 'Noh',
--[[ I don't think this is setting correctly.  cell_x is still [-1,1]^n
		mins = {0,0,0},
		maxs = {1,1,1},
--]]		
		solverVars = {
			heatCapacityRatio = 5/3,
		},
		getInitCondCode = function(self, solver)
			-- TODO add boundary condition of exact solution
			return [[
	rho = 1.;
	P = 1e-6;
	real r = real3_len(x);
	v = real3_real_mul(x, -1./r);
]]
		end,
	},
	
	-- http://www-troja.fjfi.cvut.cz/~liska/CompareEuler/compare8/node36_ct.html
	{
		name = 'implosion',
		mins = {0,0,0},
		maxs = {.3, .3, .3},
		solverVars = {
			heatCapacityRatio = 1.4,
		},
		getInitCondCode = function(self, solver)
			solver:setBoundaryMethods'mirror'
			return [[
	if (x.y < -1. - x.x) {
		rho = .125;
		P = .14;
	} else {
		rho = 1.;
		P = 1.;
	}
]]
		end,
	},
	
	-- http://www.astro.uni-bonn.de/~jmackey/jmac/node7.html
	-- http://www.astro.princeton.edu/~jstone/Athena/tests/brio-wu/Brio-Wu.html
	RiemannProblem{
		name = 'Brio-Wu',
		-- left
		{rho=1, P=1, ['B.x']=.75, ['B.y']=1, ['B.z']=0},
		-- right
		{rho=.125, P=.1, ['B.x']=.75, ['B.y']=-1, ['B.z']=0},
		solverVars = {
			heatCapacityRatio = 2,
		},
	},

	-- 2014 Abgrall, Kumar "Robust Finite Volume Scheme for Two-Fluid Plasma Equations"
	-- TODO instead of providing this as an ideal MHD initial state
	--  instead provide it as a two-fluid initial state
	-- with distinct ion and electron values
	RiemannProblem{
		name = 'two-fluid emhd modified Brio-Wu',
		{rho=1, P=5e-5, ['B.x']=.75, ['B.y']=1, ['B.z']=0},
		{rho=.125, P=5e-6, ['B.x']=.75, ['B.y']=-1, ['B.z']=0},
		solverVars = {
			heatCapacityRatio = 2,
		},
	},

	-- http://www.astro.virginia.edu/VITA/ATHENA/ot.html
	-- http://www.astro.princeton.edu/~jstone/Athena/tests/orszag-tang/pagesource.html
	-- https://www.csun.edu/~jb715473/examples/mhd2d.htm
	{
		name = 'Orszag-Tang',
		solverVars = {
			heatCapacityRatio = 5/3,
		},
		getInitCondCode = function(self, solver)
			solver:setBoundaryMethods'periodic'
			return [[
	const real B0 = 1./sqrt(4. * M_PI);
	//rho = 25./(36.*M_PI);											//Athena i.c.
	rho = solver->heatCapacityRatio * solver->heatCapacityRatio;	//CSUN i.c.
	v.x = -sin(2. * M_PI * (x.y * .5 + .5));
	v.y = sin(2. * M_PI * (x.x * .5 + .5));
	v.z = 0;
	//P = 5./(12.*M_PI);			//Athena i.c. - is this hydro pressure or total pressure?
	P = solver->heatCapacityRatio;	//CSUN i.c. - fluid pressure
	B.x = -B0 * sin(2. * M_PI * (x.y * .5 + .5));
	B.y = B0 * sin(2. * M_PI * (x.x + .5));
	B.z = 0;
]]
		end,
	},
	
	-- http://plutocode.ph.unito.it/Doxygen/Test_Problems/_m_h_d_2_rotor_2init_8c.html 
	-- http://flash.uchicago.edu/~jbgallag/2012/flash4_ug/node34.html#SECTION08123000000000000000 
	-- https://arxiv.org/pdf/0908.4362.pdf lists the following:
	--	omega = 2, P = 1, gamma = 7/5
	-- 	omega = 1, P = .5, gamma = 5/3
	{
		name = 'MHD rotor',
		solverVars = {
			heatCapacityRatio = 7/5,
		},
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			return [[
	real3 xc = coordMap(x);
	const real r0 = .1;
	const real r1 = .115;
	const real omega = 2.; 
	real r = sqrt(xc.x * xc.x + xc.y * xc.y);
	real vPhi = 0.;
	if (r <= r0) {
		rho = 10.;
		v = cartesianFromCoord(_real3(
			-omega * xc.y / r0,
			omega * xc.x / r0,
			0.
		), x);
	} else if (r <= r1) {
		real f = (r1 - r) / (r1 - r0);
		rho = 1 + 9 * f;
		v = cartesianFromCoord(_real3(
			-f * omega * xc.y / r,
			f * omega * xc.x / r,
			0.
		), x);
	} else {
		rho = 1.;
	}
	P = 1.;
	B.x = 5. / sqrt(4. * M_PI);
]]
		end,
	},

	{
		name = 'spinning magnetic fluid',
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			return [[
	real3 xc = coordMap(x);
	rho = .1;
	P = 1;
	
	real3 delta = xc;
	real coord_r = real3_len(delta);
	real3 eHat_r = real3_real_mul(delta, 1. / coord_r);
	real3 eHat_theta = _real3(-eHat_r.y, eHat_r.x, 0.);
	real3 eHat_z = _real3(0., 0., 1.);
	real radius = 1.;
	real distPastRadius = coord_r - radius;
	
	real coord_R = coord_r - .5 * (solver->mins.x + solver->maxs.x);

	if (distPastRadius < 0.) {
		rho = P = 1.;
		v = real3_real_mul(eHat_theta, .1);
#if 0
		B = real3_add(v, 
			real3_add(
				real3_real_mul(eHat_r, -x.z),
				real3_real_mul(eHat_z, coord_R)
			)
		);
#else
		B = _real3(0,0,1);
#endif
	}
]]
		end,
	},

	{
		name = 'magnetic fluid',
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			solver.useGravity = true
			return [[
	real3 xc = coordMap(x);
	rho = .1;
	P = 1;
	
	real3 delta = xc;
	real dist = real3_len(_real3(delta.x, delta.y, 0.));
	real radius = 1.;
	real distPastRadius = dist - radius;
	rho = .1;
	P = 1;
	if (distPastRadius < 0.) {
		rho = 1;
		v.z = 1;
	}
]]
		end,
	},

	{
		name = '2002 Dedner peak Bx',
		solverVars = {
			heatCapacityRatio = 5/3,
		},
		mins = {-.5, -1.5, -1},
		maxs = {.5, 1.5, 1},
		getInitCondCode = function(self, solver)
			return template([[
	rho = 1.;
	v.x = 1.;
	v.y = 1.;
	real s = 128.;	//TODO
	real r = 1. + s * s * (-128. + s * s * 4096.);
	B.x = r * (x.x * x.x + x.y * x.y) * <?=clnumber(1 / math.sqrt(4 * math.pi))?>;
	B.y = 0.;
	B.z = <?=clnumber(1 / math.sqrt(4 * math.pi))?>;
	P = 6.;
]],		{
			clnumber = clnumber,
		})
		end,
	},

	{
		name = '2002 Dedner 1D Riemann',
		solverVars = {
			heatCapacityRatio = 7/5,
		},
		mins = {-.5, -.25, -1},
		maxs = {.5, .25, 1},
		getInitCondCode = function(self, solver)
			return template([[
	rho = 1.;
	B.x = <?=clnumber(5 / math.sqrt(4 * math.pi))?>;
	B.y = <?=clnumber(5 / math.sqrt(4 * math.pi))?>;
	if (x.x < 0.) {
		v.x = 10.;
		P = 20.;
	} else {
		v.x = -10.;
		P = 1.;
	}
]], 	{
			clnumber = clnumber,
		})
		end,
	},

	{
		name = '2002 Dedner Shock Reflection',
		solverVars = {
			heatCapacityRatio = 7/5,
		},
		mins = {-1, -.5, -1},
		maxs = {1, .5, 1},
		getInitCondCode = function(self, solver)

			local boundaryMethods = {}
			for i,x in ipairs(xNames) do
				for _,minmax in ipairs(minmxas) do
					boundaryMethods[x..minmax] = 'freeflow'
				end
			end
			boundaryMethods.xmax = 'freeflow'
			boundaryMethods.ymax = 'mirror'
			solver:setBoundaryMethods(boundaryMethods)

			-- left boundary
			boundaryMethods.xmin = function(args)
				local U = 'buf['..args.index'j'..']'
				return template([[
	//TODO 'i'
	//args.index provides this ... post-converted to an integer
	//I need the vector now
	real3 x = cell_x((int4)(<?=args.indexv'j'?>,0));
	<?=eqn.prim_t?> W = {
		.rho = 1.,
		.v = _real3(2.9, 0., 0.),
		.B = _real3(.5, 0., 0.),
		.P = 5. / 7.,
<? if eqn.primStruct.vars:find(nil, function(var)
	return next(var) == 'psi'
end) then
?>		.psi = 0.,
<? end	
?>	};
	<?=U?> = consFromPrim(W, x);
]], {U=U, eqn=solver.eqn, args=args})
			end

			-- left boundary
			boundaryMethods.ymin = function(args)
				local U = 'buf['..args.index'j'..']'
				return template([[
	real3 x = cell_x(i);
	<?=eqn.prim_t?> W = {
		.rho = 1.4598,
		.v = _real3(2.717, -.4049, 0.),
		.B = _real3(.6838, -.1019, 0.),
		.P = 1.2229,
<? if eqn.primStruct.vars:find(nil, function(var)
	return next(var) == 'psi'
end) then
?>		.psi = 0.,
<? end	
?>	};
	<?=U?> = consFromPrim(W, x);
]], {U=U, eqn=solver.eqn})
			end
					
			return [[
	rho = 1.;
	v.x = 2.9;
	B.x = .5;
	P = 5./7.;
]]
		end,
	},

	{
		name = '2002 Dedner 2D Riemann problem',
		solverVars = {
			heatCapacityRatio = 5/3,
		},
		getInitCondCode = function(self, solver)
			-- TODO dirichlet 
			solver:setBoundaryMethods'freeflow'
			-- boundary: [-1, 1] x [-1, 1]
			return [[
	bool xp = x.x > mids.x;
	bool yp = x.y > mids.y;
	real3 m = real3_zero;
	real eInt = 0.;
	if (yp) {
		if (xp) {	//I
			rho = .9308;
			m = _real3(1.4557, -.4633, .0575);
			B = _real3(.3501, .9830, .3050);
			eInt = 5.0838;
		} else {	//II
			rho = 1.0304;
			m = _real3(1.5774, -1.0455, -0.1016);
			B = _real3(0.3501, 0.5078, 0.1576);
			eInt = 5.7813;
		}
	} else {
		if (!xp) {	//III
			rho = 1.0000;
			m = _real3(1.7500, -1.0000, 0.0000);
			B = _real3(0.5642, 0.5078, 0.2539);
			eInt = 6.0000;
		} else {	//IV
			rho = 1.8887;
			m = _real3(0.2334, -1.7422, 0.0733);
			B = _real3(0.5642, 0.9830, 0.4915);
			eInt = 12.999;
		}
	}
	v = real3_real_mul(m, 1. / rho);
	P = (solver->heatCapacityRatio - 1.) * rho * eInt;
]]
		end,
	},

	{
		name = '2002 Dedner Kelvin-Helmholtz',
		mins = {0,-1,-1},
		maxs = {1,1,1},
		getInitCondCode = function(self, solver)
			solver:setBoundaryMethods'periodic'
			return [[
	rho = 1.;
	v.x = 5. * (tanh(20. * (x.y + .5)) - (tanh(20. * (x.y - .5)) + 1.));
	v.y = .25 * sin(2. * M_PI * x.x) * (
		exp(-100. * (x.y + .5) * (x.y + .5))
		- exp(-100. * (x.y - .5) * (x.y - .5))
	);
	B.x = 1.;
	P = 50.;
]]
		end,
	},

	-- http://www.cfd-online.com/Wiki/Explosion_test_in_2-D
	{
		name = 'sphere',
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			return [[
	real3 xc = coordMap(x);
	real rSq = real3_lenSq(xc);
	const real R = .2;
	bool inside = rSq < R*R;
	rho = inside ? 1 : .1;
	P = inside ? 1 : .1;
]]
		end,
	},

	{
		name = 'spiral',
		guiVars = {
			{name = 'rho', value = 1},
			{name = 'P', value = 1},
			{name = 'v', value = .5},
			{name = 'D', value = 1},
		},
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			return [[
	real3 xc = coordMap(x);
	real r2 = sqrt(xc.x * xc.x + xc.y * xc.y);
	P = initCond->P;
	rho = initCond->rho;
	v.x = -xc.y * initCond->v / r2;
	v.y = xc.x * initCond->v / r2;
	D.x = -xc.y * initCond->D / r2;
	D.y = xc.x * initCond->D / r2;
]]
		end,
	},
	
	{
		name = 'cyclone',
		guiVars = {
			{name = 'rho', value = 1},
			{name = 'P', value = 1},
			{name = 'inlet_v', value = .1},
		},
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			
			local ProblemBoundary = class(solver.Boundary)
			
			ProblemBoundary.name = 'cyclone problem boundary'
		
			function ProblemBoundary:getCode(args)
				local dst
				if args.minmax == 'min' then
					dst = args.index'j'
				elseif args.minmax == 'max' then
					local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
					dst = args.index(gridSizeSide..'-numGhost+j')
				end
				local lines = table()
				if args.fields then
					for _,field in ipairs(args.fields) do
						lines:insert('buf['..dst..'].'..field..' = '..field.type..'_zero;')
					end
				else
					lines:insert(template([[
	real3 x = cell_x(i);
	real3 xc = coordMap(x);
	bool inlet = false;
	if (xc.x > 0) {
		real dy = xc.y - .5;
		real dz = z - .5;
		real dyz2 = dy*dy + dz*dy;
		inlet = dyz2 < .1;
	}
	<?=eqn.cons_t?> W = {.ptr={0}}
	W.rho = initCond->rho;
	W.v = real3_zero;
	W.P = initCond->P;
	W.ePot = 0;
	if (inlet) {
		W.v = real3(-initCond->inlet_v, 0., 0.);
	}
	buf[<?=dst?>] = primFromCons(W);
]], 				{
						dst = dst,
						eqn = solver.eqn,
					}))
				end
				return lines:concat'\n'		
			end

			local oldGetBoundaryProgramArgs = solver.getBoundaryProgramArgs
			function solver:getBoundaryProgramArgs()
				for _,xi in ipairs(xNames) do
					for _,minmax in ipairs{'min', 'max'} do
						self.boundaryMethods[xi..minmax] = ProblemBoundary()
					end
				end
				local args = oldGetBoundaryProgramArgs(self)
				args = table(args)
				return args
			end

			return [[
	rho = initCond->rho;
	P = initCond->P;
]]
		end,
	},

	{
		name = 'radial gaussian',
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			return [[
	const real gaussianCenter = 6;
	const real sigma = 1;

	real3 xc = coordMap(x);
	real r = real3_len(xc);
	real delta = (r - gaussianCenter) / sigma;
	real deltaSq = delta * delta; 
	rho = .1 + exp(-deltaSq);
	P = rho;
]]
		end,
	},
	
	{
		name = 'rarefaction wave',
		getInitCondCode = function(self, solver)
			return [[
	real delta = .1;
	rho = 1;	// lhs ? .2 : .8;
	v.x = lhs ? .5 - delta : .5 + delta;
	P = 1;
]]
		end,
	},
	
	-- 2D tests described in 2002 Kurganov, Tadmor, "Solution of Two-Dimensional Riemann Problems for Gas Dynamics without Riemann Problem Solvers"
	--  which says it is compared with  C. W. Schulz-Rinne, J. P. Collins, and H. M. Glaz, Numerical solution of the Riemann problem for two-dimensional gas dynamics
	-- and I can't find that paper right now
	quadrantProblem{
		name = 'configuration 1',
		{rho=1, P=1, ['v.x']=0, ['v.y']=0},
		{rho=.5197, P=.4, ['v.x']=-.7259, ['v.y']=0},
		{rho=.1072, P=.0439, ['v.x']=-.7259, ['v.y']=-1.4045},
		{rho=.2579, P=.15, ['v.x']=0, ['v.y']=-1.4045},
	},
	quadrantProblem{
		name = 'configuration 2',
		{rho=1, P=1, ['v.x']=0, ['v.y']=0},
		{rho=.5197, P=.4, ['v.x']=-.7259, ['v.y']=0},
		{rho=1, P=1, ['v.x']=-.7259, ['v.y']=-.7259},
		{rho=.5197, P=.4, ['v.x']=0, ['v.y']=-.7259},
	},
	quadrantProblem{
		name = 'configuration 3',
		{rho=1.5, P=1.5, ['v.x']=0, ['v.y']=0},
		{rho=.5323, P=.3, ['v.x']=1.206, ['v.y']=0},
		{rho=.138, P=.029, ['v.x']=1.206, ['v.y']=1.206},
		{rho=.5323, P=.3, ['v.x']=0, ['v.y']=1.206},
	},
	quadrantProblem{
		name = 'configuration 4',
		{rho=1.1, P=1.1, ['v.x']=0, ['v.y']=0},
		{rho=.5065, P=.35, ['v.x']=.8939, ['v.y']=0},
		{rho=1.1, P=1.1, ['v.x']=.8939, ['v.y']=.8939},
		{rho=.5065, P=.35, ['v.x']=0, ['v.y']=.8939},
	},
	quadrantProblem{
		name = 'configuration 5',
		{rho=1, P=1, ['v.x']=-.75, ['v.y']=-.5},
		{rho=2, P=1, ['v.x']=-.75, ['v.y']=.5},
		{rho=1, P=1, ['v.x']=.75, ['v.y']=.5},
		{rho=3, P=1, ['v.x']=.75, ['v.y']=-.5},
	},
	quadrantProblem{
		name = 'configuration 6',
		{rho=1, P=1, ['v.x']=.75, ['v.y']=-.5},
		{rho=2, P=1, ['v.x']=.75, ['v.y']=.5},
		{rho=1, P=1, ['v.x']=-.75, ['v.y']=.5},
		{rho=3, P=1, ['v.x']=-.75, ['v.y']=-.5},
	},
	--from SRHD Marti & Muller 2000
	{
		name = 'relativistic shock reflection',
		solverVars = {
			heatCapacityRatio = 4/3,
		},
		getInitCondCode = function(self, solver)
			solver.cfl = .5	-- needs a slower cfl
			return [[
	rho = 1;
	v.x = 1. - 1e-5;
	P = (solver->heatCapacityRatio - 1.) * rho * (1e-7 / sqrt(1. - v.x * v.x));
]]
		end,
	},
	{
		name = 'relativistic blast wave test problem 1',
		solverVars = {
			heatCapacityRatio = 5/3,
		},
		getInitCondCode = function(self, solver)
			solver.cfl = .5	-- needs a slower cfl
			return [[
	rho = lhs ? 10 : 1;
	P = (solver->heatCapacityRatio - 1.) * rho * (lhs ? 2 : 1e-6);
]]
		end,
	},
	{
		name = 'relativistic blast wave test problem 2',
		solverVars = {
			heatCapacityRatio = 5/3,
		},
		getInitCondCode = function(self, solver)
			solver.cfl = .5	-- needs a slower cfl
			return [[
	rho = 1;
	P = lhs ? 1000 : .01;
]]
		end,
	},
	{
		name = 'relativistic blast wave interaction',
		getInitCondCode = function(self, solver)
			solver.cfl = .5	-- needs a slower cfl
			return template([[
	
	bool wave1 = true
<? for i=0,solver.dim-1 do ?>
		&& x.s<?=i?> < .9 * solver->mins.s<?=i?> + .1 * solver->maxs.s<?=i?>
<? end ?>;
	bool wave2 = true
<? for i=0,solver.dim-1 do ?>
		&& x.s<?=i?> > .1 * solver->mins.s<?=i?> + .9 * solver->maxs.s<?=i?>
<? end ?>;
	rho = 1;
	P = wave1 ? 1000 : (wave2 ? 100 : .01);
]], 		{
				solver = solver,
			})
		end,
	},

	{
		name = 'Colella-Woodward',
		getInitCondCode = function(self, solver)
			solver:setBoundaryMethods'freeflow'
			return [[
	rho = 1;
	if (x.x < -.4) {
		P = 1000;
	} else if (x.x < .4) {
		P = .01;
	} else {
		P = 100;
	}
]]
		end,
	},
	
	-- derived from Athena Kelvin-Helmholtz I think
	{
		name = 'Kelvin-Helmholtz',
		
		createInitStruct = function(self, solver)
			InitCond.createInitStruct(self, solver)

			local moveAxis = 1
			local sliceAxis = 2
			
			-- move around the cylinder
			if require 'hydro.coord.cylinder'.is(solver.coord) then
				moveAxis = 2
				sliceAxis = 1
			end
			
			self:addGuiVars{
				-- these are compileTime right now
				-- also TODO initCond gui vars don't support compileTime yet
				{name='moveAxis', type='combo', value=moveAxis, options={'x','y','z'}, compileTime=true},
				{name='sliceAxis', type='combo', value=sliceAxis, options={'x','y','z'}, compileTime=true},
				{name='rhoInside', value=2.},
				{name='rhoOutside', value=1.},
				{name='amplitude', value=1e-2},
				-- not seeing much of a difference
				{name='noiseAmplitude', value=1e-4},
				{name='backgroundPressure', value=2.5},
				{name='frequency', value=2.},
				--{name='thickness', value=1e-7}
				{name='thickness', value=.025},
				{name='velInside', value=-.5},
				{name='velOutside', value=.5},
			}
		end,

		getInitCondCode = function(self, solver)	
			local boundaryMethods = {}
			for i,x in ipairs(xNames) do
				for _,minmax in ipairs(minmaxs) do
					boundaryMethods[x..minmax] = 'periodic'
					if require 'hydro.coord.cylinder'.is(solver.coord) 
					and i == 1
					then
						boundaryMethods[x..minmax] = 'mirror'
					end
				end
			end
			solver:setBoundaryMethods(boundaryMethods)
		
			return template([[
	real yq1 = solver->mins.<?=sliceAxis?> * .75 + solver->maxs.<?=sliceAxis?> * .25;
	real yq2 = solver->mins.<?=sliceAxis?> * .25 + solver->maxs.<?=sliceAxis?> * .75;

	real inside = (.5 + .5 * tanh((x.<?=sliceAxis?> - yq1) / initCond->thickness))
				- (.5 + .5 * tanh((x.<?=sliceAxis?> - yq2) / initCond->thickness));

	real theta = initCond->frequency * 2. * M_PI;
<?
for i=0,solver.dim-1 do 
	if xNames[i+1] ~= sliceAxis then
?>	theta *= (x.s<?=i?> - solver->mins.s<?=i?>) / (solver->maxs.s<?=i?> - solver->mins.s<?=i?>);
<? 
	end	
end ?>

	real noise = (solver->maxs.x - solver->mins.x) * initCond->amplitude;
	rho = inside * initCond->rhoInside + (1. - inside) * initCond->rhoOutside;
	//v.x = cos(theta) * noise;
#if dim == 2
	v.y = sin(theta) * noise;
#endif
#if dim == 3
	v.z = sin(theta) * noise;
#endif
	v.moveAxis += inside * initCond->velInside + (1. - inside) * initCond->velOutside;
	v = cartesianFromCoord(v, x);
	P = initCond->backgroundPressure;
	
	//U is initialized with random(), so use its values for unique random #s
<? assert(solver.eqn.numStates >= 5); ?>
	rho += initCond->noiseAmplitude * 2. * (U->ptr[0] - .5);
	v.x += initCond->noiseAmplitude * 2. * (U->ptr[1] - .5);
	v.y += initCond->noiseAmplitude * 2. * (U->ptr[2] - .5);
	v.z += initCond->noiseAmplitude * 2. * (U->ptr[3] - .5);
	P += initCond->noiseAmplitude * 2. * (U->ptr[4] - .5);
]],				{
					solver = solver,
					sliceAxis = self.guiVars.sliceAxis.options[self.guiVars.sliceAxis.value],
					xNames = xNames,
				}
			)
		end,
	},

	-- http://www.astro.princeton.edu/~jstone/Athena/tests/rt/rt.html
	-- TODO fixme 
	{
		name = 'Rayleigh-Taylor',
		mins = {-3,-1,-1},
		maxs = {3,1,1},
		getInitCondCode = function(self, solver)	
			local xmid = (solver.mins + solver.maxs) * .5
			
			-- triple the length along the interface dimension
			--local k = solver.dim-1
			--solver.mins.s[k] = xmid.s[k] + (solver.mins.s[k] - xmid.s[k]) * 3
			--solver.maxs.s[k] = xmid.s[k] + (solver.maxs.s[k] - xmid.s[k]) * 3

			-- triple resolution along interface dimension
			--size.s[k] = size.s[k] * 3
			-- (can it handle npo2 sizes?)

			local boundaryMethods = {}
			for i,x in ipairs(xNames) do
				for _,minmax in ipairs(minmaxs) do
					boundaryMethods[x..minmax] = i == solver.dim and 'mirror' or 'periodic'
				end
			end
			solver:setBoundaryMethods(boundaryMethods)
		
			return template([[
	const real3 externalForce = _real3(0,1,0);
	ePot = 0. <? 
for side=0,solver.dim-1 do
?> + (x.s<?=side?> - solver->mins.s<?=side?>) * externalForce.s<?=side?><?
end ?>;
	int topdim = <?=solver.dim-1?>;
	bool top = x.s[topdim] > mids.s[topdim];
	//TODO noise too?
	rho = top ? 2 : 1;
	P = 2.5 - rho * ePot;
	// or maybe it is ... pressure = (gamma - 1) * density * (2.5 - potentialEnergy)
]], {
			solver = solver,
		})
		end,
	},


	--http://www.astro.virginia.edu/VITA/ATHENA/dmr.html
	{
		name = 'double mach reflection',
		mins = {0,0,0},
		maxs = {4,1,1},
		solverVars = {
			heatCapacityRatio = 7/5,
		},
		getInitCondCode = function(self, solver)
			-- I am not correctly modeling the top boundary
			solver:setBoundaryMethods{
				xmin = 'freeflow',
				xmax = 'freeflow',
				ymin = 'mirror',
				ymax = 'freeflow',
				zmin = 'mirror',
				zmax = 'mirror',
			}
			return table{
	'#define sqrt1_3 '..clnumber(math.sqrt(1/3)),
	[[
	bool inside = x.x < x.y * sqrt1_3;
	if (inside) {
		rho = 8;
		P = 116.5;
		v.x = 8.25 * cos(30. * M_PI / 180.);
		v.y = -8.25 * sin(30. * M_PI / 180.);
	} else {
		rho = 1.4;
		P = 1;
	}
]],
			}:concat'\n'
		end,
	},

	-- http://www.cfd-online.com/Wiki/2-D_laminar/turbulent_driven_square_cavity_flow
	{
		name = 'square cavity',
		getInitCondCode = function(self, solver)
			solver:setBoundaryMethods{
				xmin = 'mirror',
				xmax = 'mirror',
				ymin = 'mirror',
				ymax = 'freeflow',	-- TODO none
				zmin = 'mirror',
				zmax = 'mirror',
			}
			return [[
	rho = 1;
	v.x = x.y > .45 ? 1 : 0;
	P = 1;
]]
		end,
	},

	{
		name = 'shock bubble interaction',
		getInitCondCode = function(self, solver)
			solver:setBoundaryMethods'freeflow'
			return [[
	const real waveX = -.45;
	real3 bubbleCenter = real3_zero;
	real bubbleRadius = .2;
	real3 delta = real3_sub(x, bubbleCenter);
	real bubbleRSq = real3_lenSq(delta);
	rho = x.x < waveX ? 1. : (bubbleRSq < bubbleRadius*bubbleRadius ? .1 : 1);
	P = x.x < waveX ? 1 : .1;
	v.x = x.x < waveX ? 0 : -.5;
]]
		end,
	},

	{
		name = 'Richmyer-Meshkov',
		mins = {-2,0,0},
		maxs = {6,1,1},
		getInitCondCode = function(self, solver)
			solver:setBoundaryMethods'freeflow'
			
			local theta = math.pi / 4
			local eta = 3
			local betaInv = .5	-- .5 for magnetic field
			local P0 = 1

			return template([[
	constant real P0 = <?=clnumber(P0)?>;
	if (x.x < -.2) {
		//TODO how do you calculate conditions for a particular mach speed shock wave? 
		//constant real gamma = heatCapacityRatio;		
		//constant real M = 2;
		//rho = (gamma + 1) / (gamma - 1) * M*M / (M*M + 2);
		//P = P0 * (2 * gamma * M*M * (gamma - 1)) / ((gamma + 1) * (gamma + 1));
		rho = 1;
		P = 10;
	} else {
		P = P0;
		if (x.x < x.y * <?=clnumber(math.tan(theta))?>) {
			rho = 1;
		} else {
			rho = <?=clnumber(eta)?>;
		}
	}
	B.x = sqrt(2 * P0 * <?=clnumber(betaInv)?>);
]], 			{
					clnumber = clnumber,
					theta = theta,
					eta = eta,
					betaInv = betaInv,
					P0 = P0,
				})
		end,
	},

	
	-- Derived from the Mara's initital condition code in Mara/test/tests.lua
	-- https://jzrake.github.io/ctf/test-page.html
	
	{
		name = 'Mara IsentropicPulse',
		guiVars = {
			{name = 'entropy_ref', value = 0.1},
			{name = 'mode', value = 2},
			{name = 'rho_ref', value = 1},
		},
		solverVars = {
			heatCapacityRatio = 1.4,
		},
		getInitCondCode = function(self, solver)
	         return [[
	real3 c = real3_real_mul(x, .5);
	real L = 1.0;
	real n = initCond->mode;
	real K = initCond->entropy_ref;
	real Gamma = solver->heatCapacityRatio;
	real rho_ref = initCond->rho_ref;
	real pre_ref = K * pow(rho_ref, Gamma);
	real cs_ref = sqrt(Gamma * pre_ref / rho_ref);
	real f = sin(n*M_PI*c.x/L);
	f *= f;
	rho = rho_ref * (1.0 + f);
	P = K * pow(rho, Gamma);
	real cs = sqrt(Gamma * P/rho);
	v.x = 2 / (Gamma - 1) * (cs - cs_ref);
]]
		end,
	},

	{
		name = 'Mara Explosion',
		guiVars = {
			{name = 'Bx', value = 4},
		},
		getInitCondCode = function(self, solver)
			return [[
	real r = coordMapR(x);
	real r2 = r * r;
	bool inside = r2 < 0.01;
	rho = inside ? 1.000 : 0.125;
	P = inside ? 1.0 : 0.1;
	B.x = initCond->Bx;
]]
		end,
	},

	{
		name = 'Mara KelvinHelmholtz',
		getInitCondCode = function(self, solver)
			solver:setBoundaryMethods'periodic'
			return template[[
	real3 c = real3_real_mul(x, .5);
	rho = fabs(c.y) > 0.25 ? 1.0 : 2.0;
	P = 2.5;
	//U is initialized with [0,1] random values	
<? assert(solver.eqn.numStates >= 2); ?>
	v.x = 0.02*(U->ptr[0] - 0.5) + fabs(c.y) > 0.25 ? -.5 : .5;
	v.y = 0.02*(U->ptr[1] - 0.5);
	v.z = 0.;
]]
		end,
	},

	{
		name = 'Mara SmoothKelvinHelmholtz',
		guiVars = {
			{name = 'P0', value = 2.5},
			{name = 'rho1', value = 1.0},
			{name = 'rho2', value = 2.0},
			{name = 'L', value = .025},
			{name = 'U1', value = 0.5},
			{name = 'U2', value = -0.5},
			{name = 'w0', value = 0.01},
		},
		getInitCondCode = function(self, solver)
			return [[
	real3 c = real3_add(real3_real_mul(x, .5), _real3(.5, .5, .5));
	const real rho1 = initCond->rho1;
	const real rho2 = initCond->rho2;
	const real L = initCond->L;
	const real U1 = initCond->U1;
	const real U2 = initCond->U2;
	const real w0 = initCond->w0;
	if (c.y < 0.25) {
		rho = rho1 - 0.5*(rho1-rho2)*exp( (c.y - 0.25)/L);
		v.x = U1 - 0.5*( U1 - U2 )*exp( (c.y - 0.25)/L);
	} else if (c.y < 0.5) {
		rho = rho2 + 0.5*(rho1-rho2)*exp(-(c.y - 0.25)/L);
		v.x = U2 + 0.5*( U1 - U2 )*exp(-(c.y - 0.25)/L);
	} else if (c.y < 0.75) {
		rho = rho2 + 0.5*(rho1-rho2)*exp( (c.y - 0.75)/L);
		v.x = U2 + 0.5*( U1 - U2 )*exp( (c.y - 0.75)/L);
	} else {
		rho = rho1 - 0.5*(rho1-rho2)*exp(-(c.y - 0.75)/L);
		v.x = U1 - 0.5*( U1 - U2 )*exp(-(c.y - 0.75)/L);
	}
	v.y = w0*sin(4.*M_PI*c.x);
	P = initCond->P0;
]]
		end,
	},

	-- TODO
	-- Mara SmoothKelvinHelmholtz
	-- Mara DensityWave
	-- Mara CollidingShocks
	-- Mara TangentialVelocity
	-- then was the Srhd cases .. I think those are Marti & Muller tests
	-- then ...
	RiemannProblem{
		name = 'Mara Shocktube1',
		{rho=1, P=1},
		{rho=.125, P=.1}},
	RiemannProblem{
		name = 'Mara Shocktube2',
		{rho=1, P=.4, ['v.x'] = -2},
		{rho=1, P=.4, ['v.x'] = 2}},
	RiemannProblem{
		name = 'Mara Shocktube3',
		{rho=1, P=1e+3},
		{rho=1, P=1e-2}},
	RiemannProblem{
		name = 'Mara Shocktube4',
		{rho=1, P=1e-2},
		{rho=1, P=1e+2}},
	RiemannProblem{
		name = 'Mara Shocktube5',
		{rho=5.99924, P=460.894, ['v.x']=19.59750},
		{rho=5.99924, P=46.095, ['v.x']=-6.19633}},
	RiemannProblem{
		name = 'Mara ContactWave',
		{rho=1.0, P=1.0, ['v.x']=0.0, ['v.y']=0.7, ['v.z']=0.2},
		{rho=0.1, P=1.0, ['v.x']=0.0, ['v.y']=0.7, ['v.z']=0.2}},
	RiemannProblem{
		name = 'Mara RMHDShocktube1',
		{rho=1.000, P=1.000, ['v.x']=0.000, ['v.y']=0.0, ['v.z']=0.0, ['B.x']=0.5, ['B.y']=1.0, ['B.z']=0.0},
		{rho=0.125, P=0.100, ['v.x']=0.000, ['v.y']=0.0, ['v.z']=0.0, ['B.x']=0.5, ['B.y']=-1.0, ['B.z']= 0.0 }}, 
	RiemannProblem{
		name = 'Mara RMHDShocktube2',
		{rho=1.080, P=0.950, ['v.x']=0.400, ['v.y']=0.3, ['v.z']=0.2, ['B.x']=2.0, ['B.y']=0.3, ['B.z']=0.3},
		{rho=1.000, P=1.000, ['v.x']=-0.450, ['v.y']=-0.2, ['v.z']=0.2, ['B.x']=2.5, ['B.y']=-0.7, ['B.z']= 0.5 }}, 
	RiemannProblem{
		name = 'Mara RMHDShocktube3',
		{rho=1.000, P=0.100, ['v.x']=0.999, ['v.y']=0.0, ['v.z']=0.0, ['B.x']=10.0, ['B.y']=0.7, ['B.z']=0.7},
		{rho=1.000, P=0.100, ['v.x']=-0.999, ['v.y']=0.0, ['v.z']=0.0, ['B.x']=10.0, ['B.y']=-0.7, ['B.z']=-0.7 }}, 
	RiemannProblem{
		name = 'Mara RMHDShocktube4',
		{rho=1.000, P=5.000, ['v.x']=0.000, ['v.y']=0.3, ['v.z']=0.4, ['B.x']=1.0, ['B.y']=6.0, ['B.z']=2.0},
		{rho=0.900, P=5.300, ['v.x']=0.000, ['v.y']=0.0, ['v.z']=0.0, ['B.x']=1.0, ['B.y']=5.0, ['B.z']= 2.0 }}, 
	RiemannProblem{
		name = 'Mara RMHDContactWave',
		{rho=10.0, P=1.0, ['v.x']=0.0, ['v.y']=0.7, ['v.z']=0.2, ['B.x']=5.0, ['B.y']=1.0, ['B.z']=0.5},
		{rho=1.0, P=1.0, ['v.x']=0.0, ['v.y']=0.7, ['v.z']=0.2, ['B.x']=5.0, ['B.y']=1.0, ['B.z']= 0.5 }}, 
	RiemannProblem{
		name = 'Mara RMHDRotationalWave',
		{rho=1, P=1, ['v.x']=0.400000, ['v.y']=-0.300000, ['v.z']=0.500000, ['B.x']=2.4, ['B.y']=1.00, ['B.z']=-1.600000},
		{rho=1, P=1, ['v.x']=0.377347, ['v.y']=-0.482389, ['v.z']=0.424190, ['B.x']=2.4, ['B.y']=-0.10, ['B.z']=-2.178213 }}, 


	-- gravity potential test - equilibrium - Rayleigh-Taylor ... or is it Jeans? (still has an shock wave ... need to fix initial conditions?)

	(function()
		local coordRadius = .5
		return {
			name = 'self-gravitation - Earth',
			-- TODO what about spherical coordinates
			mins = {-2*coordRadius, -2*coordRadius, -2*coordRadius}, 
			maxs = {2*coordRadius, 2*coordRadius, 2*coordRadius}, 
			solverVars = {
				meter = constants.EarthRadius_in_m / coordRadius,	-- radius .5, grid = 2 M_Earth, so the sphere is M_Earth size
				
				-- 4/3 pi r^3 rho = m <=> rho = 3 m / (4 pi r^3)
				kilogram = constants.EarthMass_in_kg * 3 / (4 * math.pi * coordRadius^3),
		
				-- in units of m^3/(kg s^2)
				gravitationalConstant = constants.gravitationalConstant_in_m3_per_kg_s2,
				
				coulombConstant = constants.CoulombConstant_in_kg_m3_per_C2_s2,
			},
			getInitCondCode = function(self, solver)
				local f = SelfGravProblem{
					solver = solver,
					sources={
						{
							center={0, 0, 0}, 
							radius = coordRadius,
						},
					},
				}
				return f(self, solver)
			end
		}
	end)(),

	{
		name = 'self-gravitation test 1',
		-- hmm, what about spherical coordinates...
		--mins = {-1,-1,-1},
		--maxs = {1,1,1},
		getInitCondCode = function(self, solver)
			local f = SelfGravProblem{
				solver = solver,
				sources={
					{center={0, 0, 0}, radius = .5},
				},
			}
			return f(self, solver)
		end
	},

	{
		name = 'self-gravitation test 1 spinning',
		getInitCondCode = function(self, solver)
			local inside = [[
	v.x = -2 * delta.y;
	v.y = 2 * delta.x;
	rho = 1.;
	P = 1.;
]]

			return SelfGravProblem{
				getRadiusCode = function(source)
					-- TODO compute dr's component in each dx^i's, and scale our random number by that
					-- U->rho holds noise
					return '.2 - .005 * (U->ptr[0] - .5)'
				end,
				sources={
					{
						center={0, 0, 0}, 
						radius = .2,
						-- srhd solver requires the max velocity not to exceed 1 ...
						inside = inside,
					},
				},
			}(self, solver)
		end,
	},

	{
		name = 'self-gravitation test 2',
		getInitCondCode = SelfGravProblem{ 
			sources={
				{
					center = {-.25, 0, 0},
					radius = .1,
				},
				{
					center = {.25, 0, 0},
					radius = .1,
				},
			},
		},
	},

	{
		-- TODO add tidal-locked rotations
		name = 'self-gravitation test 2 orbiting',
		getInitCondCode = SelfGravProblem{ 
			sources = {
				{
					center = {-.25, 0, 0},
					radius = .1,
					inside = [[
	v.x = -2 * delta.y - .1 * x.y;
	v.y = 2 * delta.x + .1 * x.x;
	rho = 1;
	P = 1;
					]],
				},
				{
					center = {.25, 0, 0},
					radius = .1,
					inside = [[
	v.x = -2 * delta.y - .1 * x.y;
	v.y = 2 * delta.x + .1 * x.x;
	rho = 1;
	P = 1;
					]],
				},
			},
		},
	},
	
	{
		name = 'self-gravitation test 4',
		getInitCondCode =  SelfGravProblem{
			sources={
				{center={.25, .25, 0}, radius = .1},
				{center={-.25, .25, 0}, radius = .1},
				{center={.25, -.25, 0}, radius = .1},
				{center={-.25, -.25, 0}, radius = .1},
			},
		},
	},

	{
		name = 'self-gravitation soup',
		getInitCondCode = function(self, solver)
			return [[
	int q = index;
	for (int i = 0; i < 20; ++i) {
		q = (int)floor((float)0x7fffffff * fabs(sin(2. * M_PI * (float)q / (float)0x7fffffff)));
	}
	rho = .1 * (float)(q & 0xff) / (float)0xff + .1;
	
	q = (int)floor((float)0x7fffffff * fabs(sin(2. * M_PI * (float)q / (float)0x7fffffff)));
	P = .1 * (float)(q & 0xff) / (float)0xff + .1;

	q = (int)floor((float)0x7fffffff * fabs(sin(2. * M_PI * (float)q / (float)0x7fffffff)));
	v.x = .2 * (float)(q & 0xff) / (float)0xff - .1;
	
	q = (int)floor((float)0x7fffffff * fabs(sin(2. * M_PI * (float)q / (float)0x7fffffff)));
	v.y = .2 * (float)(q & 0xff) / (float)0xff - .1;
	
	q = (int)floor((float)0x7fffffff * fabs(sin(2. * M_PI * (float)q / (float)0x7fffffff)));
	v.z = .2 * (float)(q & 0xff) / (float)0xff - .1;
]]
		end,
	},

	{
		name = 'Maxwell default',
		solverVars = {
			--meter = constants.speedOfLight_in_m_per_s,
			--speedOfLight = constants.speedOfLight_in_m_per_s,
		},
		getInitCondCode = function(self, solver)
			return template([[
	D.x = <?=eqn.susc_t?>_from_real(lhs ? 1 : -1);
	D.y = <?=eqn.susc_t?>_from_real(1);
	B.y = <?=eqn.susc_t?>_from_real(-1);
	B.z = <?=eqn.susc_t?>_from_real(lhs ? 1 : -1);
]], 		{
				eqn = solver.eqn,
			})
		end,
	},

	{
		name = 'Maxwell empty waves',
		getInitCondCode = function(self, solver)
			addMaxwellOscillatingBoundary{
				solver = solver,
				side = xNames[solver.dim]..'max',
				dir = xNames[solver.dim],
				amplitude = 1,
				period = 10,
			}
			return ''
		end,
	},

	{
		name = 'Maxwell scattering around cylinder',
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			addMaxwellOscillatingBoundary{
				solver = solver,
				side = 'xmin',
				dir = 'x',
				amplitude = 1,
				period = 10,
			}
			return template([[
	real3 xc = coordMap(x);
	if (real3_lenSq(xc) < .2*.2) {
		//2018 Balezin et al "Electromagnetic properties of the Great Pyramids..."
		permittivity = <?=eqn.susc_t?>_from_cplx(_cplx(5., .1));
	}
]], solver.eqn:getEnv())
		end,
	},

	{
		name = 'Maxwell scattering around pyramid',
		getCodePrefix = function(self, solver)
			return template([[
<? 
local table = require 'ext.table'

local objs = table()
local sqrt3 = math.sqrt(3)
if solver.dim == 2 then
	objs:insert{
		{p={0, -.5, 0}, n={0, -1, 0}},
		{p={.25*sqrt3, .25, 0}, n={.5*sqrt3, .5, 0}},
		{p={-.25*sqrt3, .25, 0}, n={-.5*sqrt3, .5, 0}},
	}
	-- if you want the floor... it is a separate CSG region
	objs:insert{
		{p={0, -.5, 0}, n={0, 1, 0}},
	}
elseif solver.dim == 3 then
	objs:insert{
		{p={0, 0, -.5}, n={0, 0, -1}},
		{p={.25*sqrt3, 0, .25}, n={.5*sqrt3, 0, .5}},
		{p={-.25*sqrt3, 0, .25}, n={-.5*sqrt3, 0, .5}},
		{p={0, .25*sqrt3, .25}, n={0, .5*sqrt3, .5}},
		{p={0, -.25*sqrt3, .25}, n={0, -.5*sqrt3, .5}},
	}
	-- if you want the floor... it is a separate CSG region
	objs:insert{
		{p={0, 0, -.5}, n={0, 0, 1}},
	}
end
?>

bool testTriangle(real3 xc) {
	return false
<? for _,obj in ipairs(objs) do ?>
		|| (true 
<? 
for _,pn in ipairs(obj) do
	local p = table(pn.p):map(clnumber):concat', '
	local n = table(pn.n):map(clnumber):concat', '
?>
			&& real3_dot(real3_sub(xc, _real3(<?=p?>)), _real3(<?=n?>)) < 0.
<? end ?>
		)
<? end ?>
	;
}
]], 		{
				solver = solver,
				clnumber = clnumber,
			})
		end,
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			-- hmm, choosing min or max doesn't matter, it always shows up on min...
			addMaxwellOscillatingBoundary{
				solver = solver,
				side = xNames[solver.dim]..'max',
				dir = xNames[solver.dim],
				amplitude = -1,
				period = 10,
			}
			return template([[
	real3 xc = coordMap(x);
	xc = real3_real_mul(xc, 2.);
	if (testTriangle(xc)) {
		//2018 Balezin et al "Electromagnetic properties of the Great Pyramids..."
		permittivity = <?=eqn.susc_t?>_from_cplx(_cplx(5., .1));
	}
]], solver.eqn:getEnv())
		end,
	},

	{
		name = 'Maxwell scattering around square',
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			-- hmm, choosing min or max doesn't matter, it always shows up on min...
			addMaxwellOscillatingBoundary{
				solver = solver,
				side = xNames[solver.dim]..'max',
				dir = xNames[solver.dim],
				amplitude = -1,
				period = 10,
			}
			return template([[
	real3 xc = coordMap(x);
	xc = real3_real_mul(xc, 2.);
	if (
		xc.x > -.5 && xc.x < .5
<? if solver.dim > 1 then ?>
		&& xc.y > -.5 && xc.y < .5
<? 	if solver.dim > 2 then ?>
		&& xc.z > -.5 && xc.z < .5
<? 	end 
end	?>
	) {
		//2018 Balezin et al "Electromagnetic properties of the Great Pyramids..."
		permittivity = <?=eqn.susc_t?>_from_cplx(_cplx(5., .1));
	}
]], solver.eqn:getEnv())
		end,
	},

	{
		name = 'Maxwell scattering around Koch snowflake',
		getCodePrefix = function(self, solver)
			return template([[

#define sqrt3 <?=clnumber(math.sqrt(3))?>

#define p1 _real3(0, .5, 0)
#define p2 _real3(.25*sqrt3, -.25, 0.)
#define p3 _real3(-.25*sqrt3, -.25, 0.)

#define n1 _real3(0, 1, 0)
#define n2 _real3(.5*sqrt3, -.5, 0)
#define n3 _real3(-.5*sqrt3, -.5, 0)

//initial branches
<? for i=1,3 do ?>
real3 branch<?=i?>(real3 x) {
	x = _real3(-x.x, -x.y, 0);	//180 rotation
	real3 n = _real3(-n<?=i?>.y, n<?=i?>.x, 0);	//angle of rotation of the normal
	x = _real3(n.x * x.x - n.y * x.y, n.x * x.y + n.y * x.x, 0.);	//rotate by 'n'
	x.y -= sqrt3*3./8.;	//translate to center
	x = real3_real_mul(x, 3.);	//scale up by 3
	x = _real3(-x.x, -x.y, 0);	//180 rotation
	return x;
}
<? end ?>

//secondary branches
real3 branch2_1(real3 x) {
	x.x += sqrt3;
	x = real3_real_mul(x, 3.);
	x.y -= .5*sqrt3;
	return x;
}

real3 branch2_2(real3 x) {
	x.x -= sqrt3;
	x = real3_real_mul(x, 3.);
	x.y -= .5*sqrt3;
	return x;
}

real3 branch2_3(real3 x) {
	real c = .5;
	real s = sqrt3*.5;
	x = _real3(x.x*c - x.y*s, x.x*s + x.y*c, 0.);	//rotate by c,s
	x = real3_real_mul(x, 3);	//scale by 3
	x.y += sqrt3;	//translate to center
	return x;
}

real3 branch2_4(real3 x) {
	real c = .5;
	real s = -sqrt3*.5;
	x = _real3(x.x*c - x.y*s, x.x*s + x.y*c, 0.);	//rotate by c,s
	x = real3_real_mul(x, 3);	//scale by 3
	x.y += sqrt3;	//translate to center
	return x;
}

bool testTriangle(real3 xc) {
	return (real3_dot(real3_sub(xc, p1), n1) < 0. &&
		real3_dot(real3_sub(xc, p2), n2) < 0. &&
		real3_dot(real3_sub(xc, p3), n3) < 0.);
}
]], 		{
				clnumber = clnumber,
			})
		end,
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			addMaxwellOscillatingBoundary{
				solver = solver,
				side = xNames[solver.dim]..'max',
				dir = xNames[solver.dim],
				amplitude = 1,
				period = 10,
			}
			
			local c = constants.speedOfLight_in_m_per_s
			local s_in_m = 1 / c
			local G = constants.gravitationalConstant_in_m3_per_kg_s2
			local kg_in_m = G / c^2
			local ke = constants.CoulombConstant_in_kg_m3_per_C2_s2
			local C_in_m = math.sqrt(ke * G) / c^2	-- m
			local Ohm_in_m = kg_in_m / (s_in_m * C_in_m^2)	-- m^0
			local resistivities = table{	-- at 20' Celsius, in Ohm m
				air = 2e+14,
				aluminum = 2.65e-8,
				copper = 1.724e-8,
				iron = 9.71e-8,
				nichrome = 1e-6,
				gold = 2.24e-8,
				silver = 1.59e-8,
				platinum = 1.06e-7,
				tungsten = 5.65e-8,
			}:map(function(v) return v * Ohm_in_m end)		
			
			return template([[
	real3 xc = coordMap(x);
	xc = real3_real_mul(xc, 2.);
	
	//conductivity = <?=clnumber(1/resistivities.air)?>;

	if (false
		|| testTriangle(xc)

<? for i=1,3 do ?>
		|| testTriangle( branch<?=i?>(xc) )	
<? end ?>

<? for i=1,3 do ?>
	<? for j=1,4 do ?>
		|| testTriangle( branch2_<?=j?>( branch<?=i?>( xc ) ) )
	<? end ?>
<? end ?>

<? for k=1,4 do ?>
	<? for j=1,4 do ?>
		<? for i=1,3 do ?>
		|| testTriangle( branch2_<?=k?>( branch2_<?=j?>( branch<?=i?>( xc ) ) ) )
		<? end ?>
	<? end ?>
<? end ?>
	
	) {
		//conductivity = 0;
		//conductivity = <?=clnumber(1/resistivities.copper)?>;
		permittivity = <?=eqn.susc_t?>_from_cplx(_cplx(5., .1));
	}

]], 		table({
				clnumber = clnumber,
				resistivities = resistivities,
			}, solver.eqn:getEnv()))
		end,
	},

	{
		name = 'Maxwell Lichtenberg',
		getInitCondCode = function(self, solver)
			local src = {math.floor(tonumber(solver.gridSize.x)*.75)}
			local dst = {math.floor(tonumber(solver.gridSize.x)*.25)}
			for j=2,solver.dim do
				src[j] = math.floor(tonumber(solver.gridSize.s[j-1])*.5)
				dst[j] = math.floor(tonumber(solver.gridSize.s[j-1])*.5)
			end
			for j=solver.dim+1,3 do
				src[j] = 0
				dst[j] = 0
			end

			local addExtraSourceProgramObj = solver.Program{
				name = 'addExtraSource',
				code = table{
					solver.modules:getCodeAndHeader(solver.sharedModulesEnabled:keys():unpack()),
					template([[
//single cell domain
kernel void addExtraSource(
	global <?=eqn.cons_t?>* UBuf
) {
	UBuf[INDEX(<?=src[1]?>,<?=src[2]?>,<?=src[3]?>)].D.x = -10;
	UBuf[INDEX(<?=dst[1]?>,<?=dst[2]?>,<?=dst[3]?>)].D.x = -10;
}
]], {
	eqn = solver.eqn,
	src = src,
	dst = dst,
}),
				}:concat'\n',
			}
			addExtraSourceProgramObj:compile()
			local addExtraSourceKernelObj = addExtraSourceProgramObj:kernel{
				name = 'addExtraSource', 
				argsOut = {solver.UBufObj},
				domain = solver.app.env:domain{dim=1, size=1},
			}
			
			local oldStep = solver.step
			function solver:step(dt)
				oldStep(self, dt)
				
				-- I just want to add to the D field at a specific point ...
				-- should this be a cpu/gpu mem cpy and write?
				-- or should this be a kernel with a single cell domain?
				addExtraSourceKernelObj()
			end
		end,
		resetState = function(self, solver)
			-- super calls applyInitCondKernel ...
			InitCond.resetState(self, solver)
			-- and here I'm going to fill the permittivity 'eps' with random noise
			-- ... and put a source + and - current 'sigma' at two points on the image
			local ptr = ffi.cast(solver.eqn.cons_t..'*', solver.UBufObj:toCPU())
			for i=0,solver.numCells-1 do
				ptr[i].sigma = math.random() * 1e-4 + 1e-7
			end
			solver.UBufObj:fromCPU(ffi.cast('real*', ptr))
		end,
	},

	{
		name = 'Maxwell wire',
		getInitCondCode = function(self, solver)
			addMaxwellOscillatingBoundary{solver=solver}
			
			local c = constants.speedOfLight_in_m_per_s
			local s_in_m = 1 / c
			local G = constants.gravitationalConstant_in_m3_per_kg_s2
			local kg_in_m = G / c^2
			local ke = constants.CoulombConstant_in_kg_m3_per_C2_s2
			local C_in_m = math.sqrt(ke * G) / c^2	-- m
			local Ohm_in_m = kg_in_m / (s_in_m * C_in_m^2)	-- m^0
			local resistivities = table{	-- at 20' Celsius, in Ohm m
				air = 2e+14,
				aluminum = 2.65e-8,
				copper = 1.724e-8,
				iron = 9.71e-8,
				nichrome = 1e-6,
				gold = 2.24e-8,
				silver = 1.59e-8,
				platinum = 1.06e-7,
				tungsten = 5.65e-8,
			}:map(function(v) return v * Ohm_in_m end)
			return template([[
	D.x = <?=eqn.susc_t?>_from_real(1.);
	
	//conductivity = <?=eqn.susc_t?>_from_real(<?=clnumber(1/resistivities.air)?>);
	
	real r2 = x.y * x.y<? if solver.dim == 3 then ?> + x.z * x.z<? end ?>;	
	
	if (r2 < .1*.1) {
		//conductivity = <?=eqn.susc_t?>_from_real(<?=clnumber(1/resistivities.copper)?>);
		permittivity = <?=eqn.susc_t?>_from_cplx(_cplx(5., .1));
#error TODO assign rhoCharge and J
	}
]], 		table({
				solver = solver,
				clnumber = clnumber,
				resistivities = resistivities,
			}, solver.eqn:getEnv()))
		end,
	},

	-- from 2000 Munz 'A Finite-Volume Method for the Maxwell Equations in the Time Domain', section 4
	{
		name = 'Maxwell transverse waves',
		guiVars = {
			{name = 'E0', value = 1, compileTime=true},
			{name = 'm', value = 8, compileTime=true},
			{name = 'n', value = 5, compileTime=true},
			{name = 'x0', value = 2, compileTime=true},
			{name = 'y0', value = 2, compileTime=true},
		},
		getInitCondCode = function(self, solver)
			return template([[
	D.z = <?=eqn.susc_t?>_from_real( E0 * sin(m * M_PI * x.x / x0) * sin(n * M_PI * x.y / y0) );
]], 		table({
				solver = solver,
			}, solver.eqn:getEnv()))
		end,
	},

	{
		guiVars = {
			{name = 'rhoCharge0', value = 1},
		},
		name = 'Maxwell charged particle',
		getInitCondCode = function(self, solver)
			return template([[
	rhoCharge = (i.x == solver->gridSize.x/2 && i.y == solver->gridSize.y/2 && i.z == solver->gridSize.z/2) ? initCond->rhoCharge0 : 0.;
]], 		table({
				solver = solver,
			}, solver.eqn:getEnv()))
		end,
	},

	{
		name = '2017 Degris et al',
		getInitCondCode = function(self, solver)
			return [[
	rho = 1.;
	P = 1.;
	if (x.x <= -.8) {
		B.x = 0.;
	} else if (x.x <= -.6) {
		B.x = -2. * (x.x + .8);
	} else if (x.x <= .6) {
		B.x = exp(-.5 * (x.x / .11) * (x.x / .11));
	} else { 
		B.x = .5;
	}
]]
		end,
	},

	{	-- from https://www.youtube.com/watch?v=Fe_f_mQCY6g
		-- doesn't have much description to it
		name = 'that one mhd simulation from youtube',
		getInitCondCode = function(self, solver)
			-- I am not correctly modeling the top boundary
			solver:setBoundaryMethods{
				xmin = 'periodic',
				xmax = 'periodic',
				ymin = 'mirror',
				ymax = 'mirror',
				zmin = 'mirror',
				zmax = 'mirror',
			}
			return [[
	rho = 1.;
	P = 1. + .1 * x.y + .1 * cos(M_PI * x.x) * cos(M_PI * .5 * x.y);
	B.x = .1;
]]
		end,
	},

	{
		name = 'spiral with flipped B field',
		guiVars = {
			{name = 'rho', value = 1},
			{name = 'P', value = 1},
			{name = 'v', value = .5},
			{name = 'B', value = 1},
		},
		depends = {'coordMap'},
		getInitCondCode = function(self, solver)
			return [[
	real3 xc = coordMap(x);
	real r = real3_len(xc);
	P = initCond->P;
	rho = initCond->rho;
	v.x = -xc.y * initCond->v / r;
	v.y = xc.x * initCond->v / r;
	real s = sign(r - .5);
	B.x = -xc.y * s * initCond->B / r;
	B.y = xc.x * s * initCond->B / r;
]]
		end,
	},

	{
		name = 'two-fluid EMHD soliton ion',
		getInitCondCode = function(self, solver)
			return [[
	const real L = 12.;
	rho = 1. + (real)exp((real)-25. * (real)fabs(x.x - L / (real)3.));
]]
		end,
	},
	{
		name = 'two-fluid EMHD soliton electron',
		getInitCondCode = function(self, solver)
			return [[
	const real L = 12.;
	rho = 5. * (1. + (real)exp((real)-25. * (real)fabs(x.x - L / (real)3.)));
]]
		end,
	},
	{
		name = 'two-fluid EMHD soliton maxwell',
		getInitCondCode = function(self, solver)
-- TODO			
			return [[

]]
		end,
	},
}:map(function(cl)
	return class(InitCond, cl)
end)
return initConds

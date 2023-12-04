--[[
initial conditions for Euler fluid equations (and similar)
slowly making this the home of all initial conditions...
--]]
local ffi = require 'ffi'
local table = require 'ext.table'
local range = require 'ext.range'
local math = require 'ext.math'
local path = require 'ext.path'
local clnumber = require 'cl.obj.number'
local symmath = require 'symmath'
local materials = require 'hydro.materials'
local InitCond = require 'hydro.init.init'
local real = require 'hydro.real'
local constants = require 'hydro.constants'

local common = require 'hydro.common'
local xNames = common.xNames
local minmaxs = common.minmaxs

local EulerInitCond = InitCond:subclass()

--[[
create mins/maxs functions that are conditional on the solver.coord
but constrained within the unit bounds of 'size'

TODO use the same techinque that the 3D volumizer uses
--]]
local function createMinsMaxs(coordRadius)
	return
		-- mins
		function(self)
			local solver = assert(self.solver)
			local coord = assert(solver.coord)
			if require 'hydro.coord.cylinder':isa(coord) then
				return {0, 0, -2*coordRadius}
			end
			if require 'hydro.coord.sphere':isa(coord) then
				return {0, 0, -math.pi}
			end
			return {-2*coordRadius, -2*coordRadius, -2*coordRadius}
		end,
		-- maxs
		function(self)
			local solver = assert(self.solver)
			local coord = assert(solver.coord)
			if require 'hydro.coord.cylinder':isa(coord) then
				return {2*coordRadius, 2*math.pi, 2*coordRadius}
			end
			if require 'hydro.coord.sphere':isa(coord) then
				return {2*coordRadius, math.pi, math.pi}
			end
			return {2*coordRadius, 2*coordRadius, 2*coordRadius}
		end
end

local function RiemannProblem(initCond)
	local WL, WR = initCond[1], initCond[2]

	local function getInitCondFieldName(i,varname)
		local suffix = ({'L','R'})[i]
		return varname:gsub('%.','')..suffix
	end

	function initCond:createInitStruct()
		EulerInitCond.createInitStruct(self)

		for i,WLR in ipairs(initCond) do
			--[[ you could just use pairs ... but then field order is arbitrary ... and then code gen can differ ... and caching is impaired ...
			for prefix, value in pairs(WLR) do
			--]]
			-- [[ so instead ....
			for _,prefix in ipairs(table.keys(WLR):sort()) do
				local value = WLR[prefix]
			--]]
				local name = getInitCondFieldName(i,prefix)
				self:addGuiVar{
					name = name,
					value =
						self.args[name] 	-- use the value specified in initCondArgs = {...},
						or value			-- use the value specified in the RiemannProblem{...} default,
				}
			end
		end
	end

	function initCond:init(args)
		EulerInitCond.init(self, args)

		if args then
			self.overrideDim = args.dim
		end
	end

	function initCond:getClassDefCode()
		local function build(i)
			return table.keys(initCond[i]):sort():mapi(function(name,_,t)
				return '\t\targs.'..name..' = args.initCond.'..getInitCondFieldName(i,name)..';', #t+1
			end):concat'\n'
		end
		return self.solver.eqn:template(
path'hydro/init/euler.clcpp':read()
..[[
namespace Hydro {
template<
	int dim_,
	typename Prim,
	typename CellArgs
>
struct InitCond_Euler_<?=initCond.name?> : public Hydro::RiemannProblem<
	dim_,
	Prim,
	CellArgs,
	InitCond_Euler_<?=initCond.name?><dim_, Prim, CellArgs>
> {
	static void buildL(CellArgs & args) {
<?=build(1)?>
	}
	static void buildR(CellArgs & args) {
<?=build(2)?>
	}
};
}

// TODO merge initCond_t i.e. InitCond with InitCondC here
namespace <?=Solver?> {
template<typename Prim, typename Cons> using InitCondC = Hydro::InitCond_Euler_<?=initCond.name?><
	<?=initCond.overrideDim or solver.dim?>,
	Prim,
	Hydro::InitCondCellArgs<Cons>
>;
}	//namespace <?=Solver?>
]], 		{
				initCond = self,
				build = build,
			})
	end


	-- http://www.itam.nsc.ru/flowlib/SRC/sod.f
	-- TODO it might be more efficient to hand this function the array, and have it return an array
	-- or at least have it calculate certain values up front before iterating across all x's
	function initCond:exactSolution(t, x, y, z)
		local solver = assert(self.solver)
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
	initCond.getInitCondCode = function(self)
		local solver = assert(self.solver)
		-- [[ specific to 2002 Kurganov, Tadmor, "Solution of Two-Dimensional Riemann Problems for Gas Dynamics without Riemann Problem Solvers"
		solver.cfl = .475
		solver:setBoundaryMethods'freeflow'
		--]]
		local function build(i)
			local q = initCond[i]
			return table.keys(q):sort():mapi(function(k)
				return k..'='..q[k]..';'
			end):concat' '
		end
		return solver.eqn:template([[
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
local SelfGravProblem = EulerInitCond:subclass()

function SelfGravProblem:getRadiusCode(source)
	return clnumber(source.radius)
end

function SelfGravProblem:init(...)
	assert(self.getRadiusCode)
	assert(self.outside)
	SelfGravProblem.super.init(self, ...)
end

function SelfGravProblem:outside()
	return [[
		rho = .1;
		P = 1;
]]
end

function SelfGravProblem:getClassDefCode()
	local solver = assert(self.solver)

	solver.useGravity = true

	return solver.eqn:template(
path'hydro/init/euler.clcpp':read()
..[[
namespace Hydro {

template<typename Prim, typename CellArgs>
struct InitCond_Euler_<?=name?> {
	static inline void initCond(CellArgs & args) {
		auto & [solver, initCond, x, U, rho, v, P, ePot, D, B] = args;
<?=self:outside()?>
		//notice about initializing random velocity -- it isn't uniform about .5 so it will pull left
		//v.x = .2 * (U.m.x - .5);	//U is initialized to random()
		//v.y = .2 * (U.m.y - .5);
		//v.z = .2 * (U.m.z - .5);

//// MODULE_DEPENDS: <?=coordMap?>
		<? for i,source in ipairs(sources) do ?>{
			real3 const xc = coordMap(x);
			real3 const delta = xc - real3(<?=
				clnumber(source.center[1])?>,<?=
				clnumber(source.center[2])?>,<?=
				clnumber(source.center[3])
			?>);
			real const distSq = lenSq(delta);
			real const radius = <?=self:getRadiusCode(source)?>;
			if (distSq < radius * radius) {
				<?=source.inside or 'rho = P = 1;'?>
			}
		}<? end ?>
	}
};
}	// namespace Hydro

namespace <?=Solver?> {
template<typename Prim, typename Cons> using InitCondC = Hydro::InitCond_Euler_<?=name?><
	Prim,
	Hydro::InitCondCellArgs<Cons>
>;
}
]], {
		self = self,
		sources = self.sources,
		clnumber = clnumber,
		name = self.name:gsub('[%- ]', '_'),
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
	local BoundaryOscillating = solver.Boundary:subclass()
	BoundaryOscillating.name = 'oscillating'
	function BoundaryOscillating:getCode(args)
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		local U = 'buf['..args.index(
			side:sub(-3) == 'min' and 'j+1' or (gridSizeSide..' - solver->numGhost - 1 - j')
		)..']'

		-- TODO put the old code here
		if not args.fields then
			return
				--oldxmin(args) ..
				solver.eqn:template([[
<?
	local epsSrc =
		(require 'hydro.eqn.twofluid-emhd-lingr':isa(eqn)
			or require 'hydro.eqn.glm-maxwell':isa(eqn)
			or require 'hydro.eqn.maxwell':isa(eqn)
		)
		and U or 'solver'
?>

	<?=U?>.B = <?=vec3?>_zero;
	<?=U?>.D = <?=vec3?>_zero;
	<?=U?>.D.<?=dir?> = <?=real_mul?>(
		<?=mul?>(<?=epsSrc?>.sqrt_1_eps, <?=epsSrc?>.sqrt_1_eps),
		<?=amplitude?> * sin(2. * M_PI * <?=frequency?> * t));
]], 		{
				U = U,
				frequency = clnumber(frequency),
				amplitude = clnumber(amplitude),
				dir = dir,
			})
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
			obj.obj:setArg(3, real(self.t))
		end
		oldBoundary(self)
	end
end


local EulerAnalytical = EulerInitCond:subclass()

function EulerAnalytical:finalizeInitStruct()
	-- createInitStruct will overwrite self.guiVars ...
	EulerAnalytical.super.finalizeInitStruct(self)
	-- ... so here is where we should add the symvars
	for _,v in ipairs(self.guiVars) do
		v.symvar = symmath.var(v.name)
	end
end

-- ok now where to do the building of the expressions?
-- how bout in getInitCondCode?
function EulerAnalytical:getInitCondCode()
	local t, x, y, z = symmath.vars('t', 'x', 'y', 'z')
	self.txVars = table{t,x,y,z}

	t:nameForExporter('C', '0')	-- only in the cl init cond, set t=0 ...
	-- assume the xyz symmath vars are supposed to be in cartesian coordinates
	for _,v in ipairs{x, y, z} do
		v:nameForExporter('C', 'xc.'..v.name)
	end
	for _,v in ipairs(self.guiVars) do
		v.symvar:nameForExporter('C', 'initCond.'..v.name)
		v.symvar:nameForExporter('Lua', 'initCondPtr.'..v.name)
	end

	self.primExprs = table{self:getPrimExprs()}
	local rhoExpr, vxExpr, vyExpr, vzExpr, PExpr = self.primExprs:unpack()

	local output = {
		{['rho'] = rhoExpr},
		{['v.x'] = vxExpr},
		{['v.y'] = vyExpr},
		{['v.z'] = vzExpr},
		{['P'] = PExpr},
	}

	local clcode = table{
		'//// MODULE_DEPENDS: <?=coordMap?>',
		'real3 const xc = coordMap(x);',
		'real const t = 0.;',
		symmath.export.C:toCode{
			assignOnly = true,
			output = output,
		}
	}:concat'\n'

	-- [[ while we're here, generate the exact solution function
	-- TODO likewise the clcode and exact solution can both be generated elsewhere
	-- and clcode just returned here

	local heatCapacityRatio = symmath.var'heatCapacityRatio'
	heatCapacityRatio:nameForExporter('Lua', 'solverPtr.heatCapacityRatio')
	local luaFuncOutput = {
		{['rho'] = rhoExpr},
		{['mx'] = rhoExpr * vxExpr},
		{['my'] = rhoExpr * vyExpr},
		{['mz'] = rhoExpr * vzExpr},
		{['ETotal'] = PExpr / (heatCapacityRatio - 1)
			+ 0.5 * rhoExpr * (vxExpr^2 + vyExpr^2 + vzExpr^2)
		},
	}

	local luafunc = symmath.export.Lua:toFunc{
		input = {symmath.var'solverPtr', symmath.var'initCondPtr', t, x, y, z},
		output = luaFuncOutput,
	}

	function self:exactSolution(t, x, y, z)
		local solver = assert(self.solver)
		local solverPtr = solver.solverPtr
		local initCondPtr = solver.initCondPtr
		return luafunc(solverPtr, initCondPtr, t, x, y, z)
	end
	--]]

	return self.solver.eqn:template(clcode)
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
		getClassDefCode = function(self)
			local args = self.args
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

			return self.solver.eqn:template(
path'hydro/init/euler.clcpp':read()
..[[
namespace <?=Solver?> {
template<typename Prim, typename Cons>
struct InitCondC {
	static inline void initCond(Hydro::InitCondCellArgs<Cons> & args) {
		auto & [solver, initCond, x, U, rho, v, P, ePot, D, B] = args;
	
		rho = initCond.rho0;
		v.x = initCond.vx0;
		v.y = initCond.vy0;
		v.z = initCond.vz0;
		P = initCond.P0;
	
	}
};
}
]])
		end,
	},
	{
		name = 'random',
		getInitCondCode = function(self)
			self.solver.useGravity = true
			return [[
	rho = U.rho + 1.;
#if dim == 2
	v = real3(
		U.m.x * cos(U.m.y * 2. * M_PI),
		U.m.x * sin(U.m.y * 2. * M_PI),
		2. * U.m.z - 1.
	);
#else
	// TODO spherical random distribution using asin and stuff
	v = real3_sub(U.m, real3(.5, .5, .5));
#endif
	P = U.ETotal + 1.;
]]
		end,
	},
	{
		name = 'linear',
		solverVars = {
			heatCapacityRatio = 7/5,
		},
		getInitCondCode = function(self)
			return self.solver.eqn:template[[
//// MODULE_DEPENDS: <?=coordMap?>
	real3 xc = coordMap(x);
	rho = 1. + xc.x;
	P = 1. + xc.x;
]]
		end,
	},

	-- 2017 Zingale "Introduction to Computational Astrophysics" section 7.9.3
	EulerAnalytical:subclass{
		name = 'advect gaussian',
		-- TODO fix the default case - it explodes
		guiVars = {
			{name = 'rho0', value = 1e-3},
			{name = 'rho1', value = 1},
			{name = 'sigma', value = .1},
			{name = 'u0', value = 0},
			{name = 'v0', value = 0},
			{name = 'w0', value = 0},
			{name = 'P0', value = 1e-6},
			{name = 'x0', value = -.5},
			{name = 'y0', value = -.5},
			{name = 'z0', value = 0},
		},
		getInitCondCode = function(self)
			self.solver:setBoundaryMethods'periodic'
			return EulerAnalytical.getInitCondCode(self)
		end,
		getPrimExprs = function(self)
			local rho0, rho1, sigma, u0, v0, w0, P0, x0, y0, z0 = self.guiVars:mapi(function(v) return v.symvar end):unpack()
			local t, x, y, z = self.txVars:unpack()

			local dx = x - x0 - u0 * t
			local dy = y - y0 - v0 * t
			local dz = z - z0 - w0 * t
			local xSq = dx^2 + dy^2 + dz^2

			local rhoExpr = (rho1 - rho0) * symmath.exp(-xSq / sigma^2) + rho0
			local vxExpr = u0
			local vyExpr = v0
			local vzExpr = w0
			local PExpr = P0

			return rhoExpr, vxExpr, vyExpr, vzExpr, PExpr
		end,
	},

	{
		-- boundary waves seem to mess with this,
		-- otherwise it looks like a wave equation solution
		name = 'Bessel',
		getInitCondCode = function(self)
			return self.solver.eqn:template[[
//// MODULE_DEPENDS: <?=coordMapR?>
	real r = coordMapR(x);
//// MODULE_DEPENDS: <?=coordMap?>
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
//// MODULE_DEPENDS: Bessel
	real vr = -BESSJ1(r);
	v = real3_real_mul(n, vr);
]]
		end,
	},

	EulerAnalytical:subclass{
		name = 'advect wave',
		mins = {0,0,0},
		maxs = {1,1,1},
		guiVars = {
			{name = 'rho0', value = 1},
			{name = 'rho1', value = 3.2e-1},
			{name = 'v0x', value = 1},
			{name = 'P0', value = 1},
		},
		solverVars = {
			heatCapacityRatio = 7/5,
		},
		getInitCondCode = function(self)
			self.solver:setBoundaryMethods'periodic'
			return EulerAnalytical.getInitCondCode(self)
		end,
		getPrimExprs = function(self)
			local rho0, rho1, v0x, P0 = self.guiVars:mapi(function(v) return v.symvar end):unpack()
			local t, x, y, z = self.txVars:unpack()

			local xmin = symmath.var'xmin'
			xmin:nameForExporter('C', 'solver->mins.x')
			xmin:nameForExporter('Lua', 'solverPtr.mins.x')

			local xmax = symmath.var'xmax'
			xmax:nameForExporter('C', 'solver->maxs.x')
			xmax:nameForExporter('Lua', 'solverPtr.maxs.x')

			local k0 = 2 * symmath.pi / (xmax - xmin)
			local rhoExpr = rho0 + rho1 * symmath.sin(k0 * (x - xmin - v0x * t))
			local vxExpr = v0x
			local vyExpr = symmath.clone(0)
			local vzExpr = symmath.clone(0)
			local PExpr = P0

			return rhoExpr, vxExpr, vyExpr, vzExpr, PExpr
		end,
	},

	-- test case vars
	-- [[ InitCond generated by RiemannProblem lua class ...
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
	--]]
	--[=[ hmm still is verbose compared to above
	{
		name = 'Sod',
		-- L = high pressure / density
		{rho = 1, P = 1},
		-- R = low pressure / density
		{rho = .125, P = .1},
		guiVars = {
			{name = 'rhoL', value = 1},
			{name = 'PL', value = 1},
			{name = 'rhoR', value = .125},
			{name = 'PR', value = .1},
		},
		solverVars = {
			heatCapacityRatio = 5/3,
		},
		getClassDefCode = function(self)
			return self.solver.eqn:template(
path'hydro/init/euler.clcpp':read()
..[[
namespace <?=Solver?> {
template<typename Prim, typename Cons> using InitCondC = Hydro::InitCond_Euler_Sod<
	<?=initCond.overrideDim or solver.dim?>,
	Prim,
	Hydro::InitCondCellArgs<Cons>
>;
}
]], 			{
					initCond = self,
				}
)
		end,
	},
	--]=]

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
		getClassDefCode = function(self)
			return self.solver.eqn:template(
path'hydro/init/euler.clcpp':read()
..[[
namespace Hydro {
template<
	typename Prim,
	typename CellArgs
>
struct InitCond_Euler_rectangle {
	static inline void initCond(
		CellArgs & args
	) {
		auto & [solver, initCond, x, U, rho, v, P, ePot, D, B] = args;

//// MODULE_DEPENDS: <?=coordMap?>
		real3 const xc = coordMap(x);
		bool const inside = true
<?
for i=1,solver.dim do
	local xi = xNames[i]
?>
			&& xc.<?=xi?> > -.5 && xc.<?=xi?> < .5
<?
end
?>
		;

		rho = inside ? initCond.rhoL : initCond.rhoR;
		P = inside ? initCond.PL : initCond.PR;
		B.x = inside ? initCond.BxL : initCond.BxR;
		B.y = inside ? initCond.ByL : initCond.ByR;
		B.z = inside ? initCond.BzL : initCond.BzR;
	}
};
}	// namespace Hydro

// TODO merge initCond_t i.e. InitCond with InitCondC here
namespace <?=Solver?> {
template<typename Prim, typename Cons> using InitCondC = Hydro::InitCond_Euler_rectangle<
	Prim,
	Hydro::InitCondCellArgs<Cons>
>;
}
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

		getInitCondCode = function(self)
			return [[
	int4 const i = globalInt4();
	rho = initCond.rho0;
	P = (i.x == solver->gridSize.x/2 && i.y == solver->gridSize.y/2 && i.z == solver->gridSize.z/2) ? initCond.P1 : initCond.P0;
]]
		end,
	},

	-- http://www-troja.fjfi.cvut.cz/~liska/CompareEuler/compare8/node36_ct.html
	{
		name = 'Noh',
--[[ I don't think this is setting correctly.  cell pos is still [-1,1]^n
		mins = {0,0,0},
		maxs = {1,1,1},
--]]
		solverVars = {
			heatCapacityRatio = 5/3,
		},
		getInitCondCode = function(self)
			-- TODO add boundary condition of exact solution
			return [[
	rho = 1.;
	P = 1e-6;
	real r = x.length();
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
		getInitCondCode = function(self)
			self.solver:setBoundaryMethods'mirror'
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
		getInitCondCode = function(self)
			self.solver:setBoundaryMethods'periodic'
			return [[
	real const B0 = 1./sqrt(4. * M_PI);
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
		getInitCondCode = function(self)
			return self.solver.eqn:template[[
//// MODULE_DEPENDS: <?=coordMap?>
	real3 xc = coordMap(x);
	real const r0 = .1;
	real const r1 = .115;
	real const omega = 2.;
	real r = sqrt(xc.x * xc.x + xc.y * xc.y);
	real vPhi = 0.;
//// MODULE_DEPENDS: <?=cartesianFromCoord?>
	if (r <= r0) {
		rho = 10.;
		v = cartesianFromCoord(real3(
			-omega * xc.y / r0,
			omega * xc.x / r0,
			0.
		), x);
	} else if (r <= r1) {
		real f = (r1 - r) / (r1 - r0);
		rho = 1 + 9 * f;
		v = cartesianFromCoord(real3(
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
		getInitCondCode = function(self)
			return self.solver.eqn:template[[
//// MODULE_DEPENDS: <?=coordMap?>
	real3 xc = coordMap(x);
	rho = .1;
	P = 1;

	real3 delta = xc;
	real coord_r = delta.length();
	real3 eHat_r = real3_real_mul(delta, 1. / coord_r);
	real3 eHat_theta = real3(-eHat_r.y, eHat_r.x, 0.);
	real3 eHat_z = real3(0., 0., 1.);
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
		B = real3(0,0,1);
#endif
	}
]]
		end,
	},

	{
		name = 'magnetic fluid',
		getInitCondCode = function(self)
			self.solver.useGravity = true
			return self.solver.eqn:template[[
//// MODULE_DEPENDS: <?=coordMap?>
	real3 xc = coordMap(x);
	rho = .1;
	P = 1;

	real3 delta = xc;
	real dist = real3(delta.x, delta.y, 0.).length();
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
		getInitCondCode = function(self)
			return self.solver.eqn:template([[
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
		getInitCondCode = function(self)
			return self.solver.eqn:template([[
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
		getInitCondCode = function(self)
			local solver = assert(self.solver)
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
				return solver.eqn:template([[
	//TODO 'i'
	//args.index provides this ... post-converted to an integer
	//I need the vector now
	real3 const x = cellBuf[<?=args.index'j'?>].pos;
	<?=prim_t?> W = {
		.rho = 1.,
		.v = real3(2.9, 0., 0.),
		.B = real3(.5, 0., 0.),
		.P = 5. / 7.,
<? if eqn.primStruct.vars:find(nil, function(var)
	return next(var) == 'psi'
end) then
?>		.psi = 0.,
<? end
?>	};
	<?=U?> = consFromPrim(W, x);
]], {U=U, args=args})
			end

			-- left boundary
			boundaryMethods.ymin = function(args)
				local U = 'buf['..args.index'j'..']'
				return solver.eqn:template([[
	real3 const x = cellBuf[index].pos;
	<?=prim_t?> W = {
		.rho = 1.4598,
		.v = real3(2.717, -.4049, 0.),
		.B = real3(.6838, -.1019, 0.),
		.P = 1.2229,
<? if eqn.primStruct.vars:find(nil, function(var)
	return next(var) == 'psi'
end) then
?>		.psi = 0.,
<? end
?>	};
	<?=U?> = consFromPrim(W, x);
]], {U=U})
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
		getInitCondCode = function(self)
			-- TODO dirichlet
			self.solver:setBoundaryMethods'freeflow'
			-- boundary: [-1, 1] x [-1, 1]
			return [[
	bool xp = x.x > mids.x;
	bool yp = x.y > mids.y;
	real3 m = real3_zero;
	real eInt = 0.;
	if (yp) {
		if (xp) {	//I
			rho = .9308;
			m = real3(1.4557, -.4633, .0575);
			B = real3(.3501, .9830, .3050);
			eInt = 5.0838;
		} else {	//II
			rho = 1.0304;
			m = real3(1.5774, -1.0455, -0.1016);
			B = real3(0.3501, 0.5078, 0.1576);
			eInt = 5.7813;
		}
	} else {
		if (!xp) {	//III
			rho = 1.0000;
			m = real3(1.7500, -1.0000, 0.0000);
			B = real3(0.5642, 0.5078, 0.2539);
			eInt = 6.0000;
		} else {	//IV
			rho = 1.8887;
			m = real3(0.2334, -1.7422, 0.0733);
			B = real3(0.5642, 0.9830, 0.4915);
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
		getInitCondCode = function(self)
			self.solver:setBoundaryMethods'periodic'
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
		createInitStruct = function(self)
			EulerInitCond.createInitStruct(self)
			local args = self.args
			self:addGuiVars{
				{name = 'radius', value = args.radius or .5},
				{name = 'rhoInside', value = args.rhoInside or 1},
				{name = 'PInside', value = args.rhoInside or 1},
				{name = 'rhoOutside', value = args.rhoOutside or .01},
				{name = 'POutside', value = args.rhoOutside or .01},
			}
		end,
		getInitCondCode = function(self)
			return self.solver.eqn:template[[
//// MODULE_DEPENDS: <?=coordMap?>
	real3 const xc = coordMap(x);
	real const rSq = real3_lenSq(xc);
	real const radius = initCond.radius;
	bool const inside = rSq < radius*radius;
	rho = inside ? initCond.rhoInside : initCond.rhoOutside;
	P = inside ? initCond.PInside : initCond.POutside;
]]
		end,
	},

	{
		name = 'spiral',
		createInitStruct = function(self)
			EulerInitCond.createInitStruct(self)
			local args = self.args
			-- TODO instead of this here, it'd be nice to have a static tuple fields
			-- and then ... hmm how would the gui iterate over it ... there still has to be some lua reflection for the gui to work
			-- even if I wrote the non-CL code in C++ , it' still need some CL<->C++ reflection layer code...
			-- so it is easiest to put it in Lua here (or C++ ofc ... if only C++ had as easy of string and table manipulation as Lua does)
			self:addGuiVars{
				{name = 'torusGreaterRadius', value = args.torusGreaterRadius or .5},
				{name = 'torusLesserRadius', value = args.torusLesserRadius or .1},	-- unitless, i.e. coordinate units
				{name = 'rhoInside', value = 1},
				{name = 'rhoOutside', value = 1e-3},
				{name = 'P', value = 1},
				{name = 'v', value = args.v or .5, units = 'm'},
				{name = 'D', value = 1, units = 'C/m^2'},
			}
		end,
		getClassDefCode = function(self)
			return self.solver.eqn:template(
path'hydro/init/euler.clcpp':read()
..[[
namespace <?=Solver?> {
template<typename Prim, typename Cons> using InitCondC = InitCond_Euler_spiral<
	<?=overrideDim or solver.dim?>,
	Prim,
	Hydro::InitCondCellArgs<Cons>
>; 
}
]],			{
				overrideDim = self.overrideDim,
			})
		end,
	},

	{
		name = 'jet',
		init = function(self, args)
			InitCond.init(self, args)

			-- add to solver, not self, so boundary can read them
			-- TODO should we pass initCond to boundary?
			self.solver.eqn:addGuiVars{
				{name = 'init_rho', value = .1},
				{name = 'init_P', value = 1},
				{name = 'init_inlet_rho', value = 1},
				{name = 'init_inlet_v', value = .1},
				{name = 'init_inlet_P', value = 1},
				{name = 'inlet_r', value = .02},
			}
		end,
		getInitCondCode = function(self)
			local solver = assert(self.solver)

			local ProblemBoundary = solver.Boundary:subclass()

			ProblemBoundary.name = 'jet boundary'

			function ProblemBoundary:getCode(args)
				local dst
				if args.minmax == 'min' then
					dst = args.index'j'
				elseif args.minmax == 'max' then
					local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
					dst = args.index(gridSizeSide..' - solver->numGhost + j')
				end
				local lines = table()
				if args.fields then
					for _,field in ipairs(args.fields) do
						lines:insert('buf['..dst..'].'..field..' = 0.;')
					end
				else
					lines:insert(solver.eqn:template([[
<? local isSRHD = require 'hydro.eqn.srhd':isa(eqn) ?>
{
	real3 const x = cellBuf[<?=dst?>].pos;
//// MODULE_DEPENDS: <?=coordMap?>
	real3 const xc = coordMap(x);
	bool inlet = false;
	if (xc.x < 0) {
		real dy = xc.y;// - .5;
		real dz = <?= solver.dim == 3 and 'xc.z - .5' or '0'?>;
		real dyzSq = dy*dy + dz*dz;
		inlet = dyzSq < solver->inlet_r * solver->inlet_r;
	}
	<?=isSRHD and prim_only_t or prim_t?> prim = {.s={0}};

	if (inlet) {
		prim.rho = solver->init_inlet_rho;
		prim.v = real3(solver->init_inlet_v, 0., 0.);
<? if not isSRHD then ?>
		prim.P = solver->init_inlet_P;
		prim.ePot = 0;
<? else ?>
		prim.eInt = calc_eInt_from_P(solver, prim.rho, solver->init_inlet_P);
<? end ?>


<? if isSRHD then ?>
		consFromPrimOnly(buf + <?=dst?>, solver, &prim, x);
<? else ?>
		<?=consFromPrim?>(buf + <?=dst?>, solver, &prim, x);
<? end ?>

	} else {

		if (xc.x > 0) {
//freeflow b.c.
<?
do
	local dst, src
	if args.minmax == 'min' then
		dst = args.index'j'
		src = args.index'solver->numGhost'
	elseif args.minmax == 'max' then
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		dst = args.index(gridSizeSide..' - solver->numGhost + j')
		src = args.index(gridSizeSide..' - solver->numGhost - 1')
	end
?>	<?=bc:assignDstSrc(dst, src, args)?>
<?
end
?>
		} else {

//constant b.c. around the inlet
		prim.rho = solver->init_rho;
		prim.v = real3_zero;
<? if not isSRHD then ?>
		prim.P = solver->init_P;
		prim.ePot = 0;
<? else ?>
		prim.eInt = calc_eInt_from_P(solver, prim.rho, solver->init_P);
<? end ?>

<? if isSRHD then ?>
		consFromPrimOnly(buf + <?=dst?>, solver, &prim, x);
<? else ?>
		<?=consFromPrim?>(buf + <?=dst?>, solver, &prim, x);
<? end ?>
		}
	}

}
]], 				{
						bc = self,
						args = args,
						dst = dst,
						eqn = solver.eqn,
					}))
				end
				local depends = table{
					assert(solver.eqn.symbols.prim_t),
					assert(solver.eqn.symbols.consFromPrim),
					assert(solver.coord.symbols.coordMap),
				}
				if require 'hydro.eqn.srhd':isa(solver.eqn) then
					depends:append{
						assert(solver.eqn.symbols.eqn_common),
					}
				end

				return lines:concat'\n', depends
			end

			local oldGetBoundaryProgramArgs = solver.getBoundaryProgramArgs
			function solver:getBoundaryProgramArgs()
				for _,xi in ipairs(xNames) do
					for _,minmax in ipairs{'min', 'max'} do
						local side = xi..minmax
						if (solver.coord.name == 'cylinder' and side ~= 'xmax')
						or (solver.coord.name == 'spherical' and side ~= 'xmax')
						then
							-- don't replace
						else
							self.boundaryMethods[side] = ProblemBoundary()
						end
					end
				end
				local args = oldGetBoundaryProgramArgs(self)
				args = table(args)
				return args
			end

			return [[
	rho = solver->init_rho;
	P = solver->init_P;
]]
		end,
	},

	{
		name = 'radial gaussian',
		getInitCondCode = function(self)
			return self.solver.eqn:template[[
	real const gaussianCenter = 6;
	real const sigma = 1;

//// MODULE_DEPENDS: <?=coordMap?>
	real3 xc = coordMap(x);
	real r = xc.length();
	real delta = (r - gaussianCenter) / sigma;
	real deltaSq = delta * delta;
	rho = .1 + exp(-deltaSq);
	P = rho;
]]
		end,
	},

	{
		name = 'rarefaction wave',
		getInitCondCode = function(self)
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
		getInitCondCode = function(self)
			self.solver.cfl = .5	-- needs a slower cfl
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
		getInitCondCode = function(self)
			self.solver.cfl = .5	-- needs a slower cfl
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
		getInitCondCode = function(self)
			self.solver.cfl = .5	-- needs a slower cfl
			return [[
	rho = 1;
	P = lhs ? 1000 : .01;
]]
		end,
	},
	{
		name = 'relativistic blast wave interaction',
		getInitCondCode = function(self)
			self.solver.cfl = .5	-- needs a slower cfl
			return self.solver.eqn:template([[

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
		getInitCondCode = function(self)
			self.solver:setBoundaryMethods'freeflow'
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
		createInitStruct = function(self)
			EulerInitCond.createInitStruct(self)
			local args = self.args or {}

			local moveAxis = 1
			local sliceAxis = 2

			-- move around the cylinder
			if require 'hydro.coord.cylinder':isa(self.solver.coord) then
				moveAxis = 2
				sliceAxis = 1
			end

			moveAxis = args.moveAxis or moveAxis
			sliceAxis = args.sliceAxis or sliceAxis

			self:addGuiVars{
				-- these are compileTime right now
				-- also TODO initCond gui vars don't support compileTime yet
				{name = 'moveAxis', type = 'combo', value = moveAxis, options = {'x','y','z'}, compileTime = true},
				{name = 'sliceAxis', type = 'combo', value = sliceAxis, options = {'x','y','z'}, compileTime = true},
				{name = 'rhoInside', value = args.rhoInside or 2.},
				{name = 'rhoOutside', value = args.rhoOutside or 1.},
				{name = 'amplitude', value = args.amplitude or 1e-2},
				-- not seeing much of a difference
				{name = 'noiseAmplitude', value = args.noiseAmplitude or 1e-2},
				{name = 'backgroundPressure', value = args.backgroundPressure or 2.5},
				{name = 'frequency', value = args.frequency or 2.},
				--{name = 'thickness', value = 1e-7}
				{name = 'thickness', value = args.thickness or .025},
				{name = 'velInside', value = args.velInside or -.5},
				{name = 'velOutside', value = args.velOutside or .5},
			}
		end,
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			local boundaryMethods = {}
			for i,x in ipairs(xNames) do
				for _,minmax in ipairs(minmaxs) do
					boundaryMethods[x..minmax] = 'periodic'
					if require 'hydro.coord.cylinder':isa(solver.coord)
					and i == 1
					then
						boundaryMethods[x..minmax] = 'mirror'
					end
				end
			end
			solver:setBoundaryMethods(boundaryMethods)

			return solver.eqn:template([[
	real yq1 = solver->mins.<?=sliceAxis?> * .75 + solver->maxs.<?=sliceAxis?> * .25;
	real yq2 = solver->mins.<?=sliceAxis?> * .25 + solver->maxs.<?=sliceAxis?> * .75;

	real inside = (.5 + .5 * tanh((x.<?=sliceAxis?> - yq1) / initCond.thickness))
				- (.5 + .5 * tanh((x.<?=sliceAxis?> - yq2) / initCond.thickness));

	real theta = initCond.frequency * 2. * M_PI;
<?
for i=0,solver.dim-1 do
	if xNames[i+1] ~= sliceAxis then
?>	theta *= (x.s<?=i?> - solver->mins.s<?=i?>) / (solver->maxs.s<?=i?> - solver->mins.s<?=i?>);
<?
	end
end ?>

#if dim == 2
#define perpAxis y
#elif dim == 3
#define perpAxis z
#endif

	real noise = (solver->maxs.x - solver->mins.x) * initCond.amplitude;
	rho = inside * initCond.rhoInside + (1. - inside) * initCond.rhoOutside;
	//v.x = cos(theta) * noise;
#if dim >= 2
	v.perpAxis = sin(theta) * noise;
#endif
	v.moveAxis += inside * initCond.velInside + (1. - inside) * initCond.velOutside;
//// MODULE_DEPENDS: <?=cartesianFromCoord?>
	v = cartesianFromCoord(v, x);
	P = initCond.backgroundPressure;

	//U is initialized with random(), so use its values for unique random #s
<? assert(solver.eqn.numStates >= 5); ?>
	rho += initCond.noiseAmplitude * 2. * (U.s[0] - .5);
#if 0
	v.x += initCond.noiseAmplitude * 2. * (U.s[1] - .5);
	v.y += initCond.noiseAmplitude * 2. * (U.s[2] - .5);
	v.z += initCond.noiseAmplitude * 2. * (U.s[3] - .5);
#elif dim >= 2
	real noisePhi = 2. * M_PI * U.s[1];
	real noiseR = U.s[2];
	v.moveAxis += initCond.noiseAmplitude * noiseR * cos(noisePhi);
	v.perpAxis += initCond.noiseAmplitude * noiseR * sin(noisePhi);
#endif
	P += initCond.noiseAmplitude * 2. * (U.s[4] - .5);
]],				{
					sliceAxis = self.guiVars.sliceAxis:getValue(),
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
		getInitCondCode = function(self)
			local solver = assert(self.solver)
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

			return solver.eqn:template([[
	real3 const externalForce = real3(0,1,0);
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
]])
		end,
	},

	{
		name = 'Taylor-Green',
		getInitCondCode = function(self)
			return self.solver.eqn:template[[
//// MODULE_DEPENDS: <?=coordMap?>
	real3 const xc = coordMap(x);
	real const ux = xc.x * 2. * M_PI;
	real const uy = xc.y * 2. * M_PI;
	rho = 1.;
	v.x = sin(ux) * cos(uy);
	v.y = -cos(ux) * sin(uy);
	P = 100. / solver->heatCapacityRatio + .25 * (cos(2. * ux) + cos(2. * uy));
]]
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
		getInitCondCode = function(self)
			local solver = assert(self.solver)
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

	--[[
	http://www.cfd-online.com/Wiki/2-D_laminar/turbulent_driven_square_cavity_flow
	initVel=2 works, initVel=10 explodes
	crashes for incompressible (maybe due to the b.c.)
	--]]
	{
		name = 'square cavity',
		guiVars = {
			{name='initVel', value=2, compileTime=true},
		},
		getInitCondCode = function(self)
			local solver = assert(self.solver)

			solver:setBoundaryMethods{
				xmin = 'mirror',
				xmax = 'mirror',
				ymin = 'mirror',

				ymax = {
					name = 'fixed',
					args = {
						fixedCode = function(self, args, dst)
							return args.solver.eqn:template([[
<?=prim_t?> W = {
	.rho = 1.,
	.P = 1.,
	.v = real3(initVel, 0., 0.),
	.ePot = 0.,
};
<?=consFromPrim?>(buf + <?=dst?>, solver, &W, cellBuf[<?=dst?>].pos);
]], 						{
								args = args,
								dst = dst,
							}), {
								args.solver.eqn.symbols.consFromPrim,
								args.solver.eqn.symbols.initCond_guiVars_compileTime,
							}
						end,
					},
				},

				zmin = 'mirror',
				zmax = 'mirror',
			}
			return [[
	rho = 1;
	P = 1;
]]
		end,
	},

	{
		name = 'shock bubble interaction',
		guiVars = {
			{name='shockwaveAxis', value=0},
			{name='waveX', value=-.45},
			{name='bubbleRadius', value=.2},
			{name='bubbleCenterX', value=0},
			{name='bubbleCenterY', value=0},
			{name='bubbleCenterZ', value=0},
			{name='rhoL', value=1},			-- before the shock wave
			{name='PL', value=1},
			{name='vL', value=0},
			{name='rhoR', value=1},			-- after the shock wave
			{name='PR', value=.1},
			{name='vR', value=-.5},
			{name='rhoInside', value=.1},	-- inside the bubble
			{name='PInside', value=.1},
			{name='vInside', value=-.5},
		},
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			solver:setBoundaryMethods'freeflow'
			return self.solver.eqn:template[[
//// MODULE_DEPENDS: <?=coordMap?>
real3 const xc = coordMap(x);
real3 const bubbleCenter = real3(initCond.bubbleCenterX, initCond.bubbleCenterY, initCond.bubbleCenterZ);
real const bubbleRadiusSq = initCond.bubbleRadius * initCond.bubbleRadius;
real3 const delta = real3_sub(xc, bubbleCenter);
real const bubbleRSq = real3_lenSq(delta);
int const axis = initCond.shockwaveAxis;
if (bubbleRSq < bubbleRadiusSq) {
	rho = initCond.rhoInside;
	P = initCond.PInside;
	v.s[axis] = initCond.vInside;
} else if (xc.s[axis] < initCond.waveX) {
	rho = initCond.rhoL;
	P = initCond.PL;
	v.s[axis] = initCond.vL;
} else {
	rho = initCond.rhoR;
	P = initCond.PR;
	v.s[axis] = initCond.vR;
}
]]
		end,
	},

	{
		name = 'Richmyer-Meshkov',
		mins = {-2,0,0},
		maxs = {6,1,1},
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			solver:setBoundaryMethods'freeflow'

			local theta = math.pi / 4
			local eta = 3
			local betaInv = .5	-- .5 for magnetic field
			local P0 = 1

			return solver.eqn:template([[
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
		getInitCondCode = function(self)
			local solver = assert(self.solver)
	         return [[
	real3 c = real3_real_mul(x, .5);
	real L = 1.0;
	real n = initCond.mode;
	real K = initCond.entropy_ref;
	real Gamma = solver->heatCapacityRatio;
	real rho_ref = initCond.rho_ref;
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
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			return [[
	real r = coordMapR(x);
	real r2 = r * r;
	bool inside = r2 < 0.01;
	rho = inside ? 1.000 : 0.125;
	P = inside ? 1.0 : 0.1;
	B.x = initCond.Bx;
]]
		end,
	},

	{
		name = 'Mara KelvinHelmholtz',
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			solver:setBoundaryMethods'periodic'
			return solver.eqn:template[[
	real3 c = real3_real_mul(x, .5);
	rho = fabs(c.y) > 0.25 ? 1.0 : 2.0;
	P = 2.5;
	//U is initialized with [0,1] random values
<? assert(solver.eqn.numStates >= 2); ?>
	v.x = 0.02*(U.s[0] - 0.5) + fabs(c.y) > 0.25 ? -.5 : .5;
	v.y = 0.02*(U.s[1] - 0.5);
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
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			return [[
	real3 c = real3_add(real3_real_mul(x, .5), real3(.5, .5, .5));
	real const rho1 = initCond.rho1;
	real const rho2 = initCond.rho2;
	real const L = initCond.L;
	real const U1 = initCond.U1;
	real const U2 = initCond.U2;
	real const w0 = initCond.w0;
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
	P = initCond.P0;
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

		-- radius .5, grid = 2 M_Earth, so the sphere is M_Earth size
		local meter = constants.EarthRadius_in_m / coordRadius

		-- With this, setting rho=1 will cause the total mass of the grid to be the mass of the Earth.
		-- 4/3 pi r^3 rho = m <=> rho = 3 m / (4 pi r^3)
		local sphereCoordVolume = (4/3) * math.pi * coordRadius * coordRadius * coordRadius
		local kilogram = constants.EarthMass_in_kg / sphereCoordVolume

		-- keep the speed of light at 1?
		--local second = meter / constants.speedOfLight_in_m_per_s
		-- but this puts our simulation's 1 unit of time at 0.05 seconds, so the simulation runs slow.  20x slower than if we just kept second at its default.
		-- so how about speeding it up?
		-- speedOfLight / (meter / second) = 1
		-- second = meter / speedOfLight
		local second = 1		-- converges steady enough
		--local second = 60*60*24
		--local second = meter / constants.speedOfLight_in_m_per_s	-- converges slow, oscillates a bit

		local mins, maxs = createMinsMaxs(coordRadius)
		return {
			name = 'self-gravitation - Earth',

			mins = mins,
			maxs = maxs,
			solverVars = {
				meter = meter,
				kilogram = kilogram,
				second = second,

				speedOfLight = constants.speedOfLight_in_m_per_s,
				divPsiWavespeed_g = constants.speedOfLight_in_m_per_s,
				divPhiWavespeed_g = constants.speedOfLight_in_m_per_s,

				gravitationalConstant = constants.gravitationalConstant_in_m3_per_kg_s2,
				coulombConstant = constants.CoulombConstant_in_kg_m3_per_C2_s2,
			},
			getInitCondCode = function(self)
				return self.solver.eqn:template([[
//// MODULE_DEPENDS: <?=coordMap?>
	real3 const xc = coordMap(x);
//// MODULE_DEPENDS: <?=coordMapR?>
	real const r = coordMapR(x);
	real const rSq = r * r;
	real const R = <?=clnumber(coordRadius)?>;
	if (rSq < R * R) {
		rho = 1.;
	} else {
		rho = 1e-3;
	}
	P = 1.;

	real const omegaLen_in_1_day = 1.;
	// 1 rev / 1 sidereal day * 1 sidereal day / 86164.0905 second
	real const omegaLen_in_1_s = omegaLen_in_1_day / 86164.0905;
	// units of 1/s
	real3 const omega = real3(0, 0, omegaLen_in_1_s);

	// v = omega cross x (in Cartesian coordinates)
	v = real3_cross(omega, xc);
]], 			{
					coordRadius = coordRadius,
				})
			end
		}
	end)(),

	-- 2021 Ludwig "Galactic rotation curve and dark matter according to gravitomagnetism"
	--[[
	-- just before eqn D.12
	L_sun = 3.828e+26 -- [W]

	parsec_in_m = 648000 / M_PI * 149597870700
	kpc = 1e+3 * parsec_in_m

	-- section 7:
	mu_0 = 22.27
	alpha_eff = 116.5		-- effective half-angle in arcseconds
	r_eff = 1.69 			-- [kpc]
	alpha_0 = 140.3
	r_0 = 2.04 				-- [kpc] = used for calculating normalized mass density profile
							-- r_0 and alpha_0 should be related by r_0 = d (pi / (180 * 3600)) * alpha_0
	s_0 = 0.360
	b_1 = 0.00245			-- b_s ? coefficients of Sersec profile for flux density
	b_2 = 0.0000344
	b_3 = -1.41e-7
	b_4 = -4.05e-10
	d_2 = -3.22e-6
	d_3 = -1.11e-9
	d_4 = 5.73e-11
	alpha_e = 353.0
	s_e = 0.874
	M_s = -15.3			-- absolute magnitude .. in what units? unitless?
	L_s = 1.02e+8		-- [L_sun] = partial radiatn flux
	m_d = 12.1			-- = apparent magnitude ... in what units?
	r_max = 12.2 		-- [kpc] = maximum galactic radius
	l = 3
	r_s = 1.46e-6		-- [kpc]
	a = 7.19 			-- [kpc]
	b = 0.567 			-- [kpc]
	r_max = 12.2 		-- [kpc]
	l_beta = 8.29		-- [kpc]
	M = 1.52e+10		-- [M_sun]
	lambda = 0.134
	rho_0 = 3.31e-20	-- [kg/m^3]
	Upsilon = 150 		-- [Upsilon_sun] = total mass-to-light ratio

	-- eqn D.11
	mu(alpha) = mu_0 +
		0 <= alpha and alpha <= alpha_0 and
			5 / (2 * log(10)) * pow(alpha / alpha_1, 1 / s_1)
		or
			5 / (2 * log(10)) * pow(alpha_0 / alpha_1, 1 / s_1) * (1 - s_2 / s_1 + s_2 / s_1 * pow(alpha / alpha_0, 1 / s_2))

	-- eqn D.19
	s(alpha) = log(alpha / alpha_eff) / log(2 * log(10) / 5 * (mu(alpha) - mu_0))
	-- Does this mean that s_0 = s(0) ? and s_1 = s(1), s_2 = s(2) ?

	-- What is alpha_1 ?
	-- After D.11, the alpha_1 for NGC 1560 is defined as ...
	alpha_1 = 99.05 -- [arcsec]
	-- ... but ... this is in the appendix.  NGC 1560 is defined in section 7.
	-- So what is the alpha_1 of the other galaxies in the paper?
	-- Also in the appendix D paragraph on variables for NGC 1560, it lists ...
	-- mu_0 = 22.28, alpha_0 = 61.46 arcsec, s_1 = 0.435, s_2 = 1.144.
	-- Why don't these match section 7's variables on NGC 1560?

	-- just after eqn D.12:
	r_1 = d * (M_PI / (180 * 3600)) * alpha_1

	-- eqn D.12
	rhoNormalized(r, 0) = 10^(-.4 * (mu(r) - mu_0))
	-- or defined as ...
	if r <= r_0 then
		rhoNormalized = exp(-pow(r/r_1, 1/s_1))
	else
		rhoNormalized = exp(-pow(r_0/r_1, 1/s_1) * (1 - s_2/s_1 + s_2/s_1 * pow(r/r_0, 1/s_2)))
	end
	-- What units is this equation?
	-- What is s_1 and s_2?  Related to s_0 or s_e?
	-- Or defined in Appendix D's info on NGC 1560, a paragraph where other vars listed for NGC 1560 don't match section 7's vars of NGC 1560.

	-- eqn C.4
	rho = rhoNormalized * R_0 * R_0 * R_0 * rho_0

	-- what is R_0 ?
	-- from inline before 4.5 and from inline before C.3
	-- why is R_0 = 1 kpc?  esp when the flux function uses R_0 as a cutoff.
	R_0 = 1 -- [kpc]

	-- what is rho_0?  is (somewhere) rho(0,0) at the center of the function.

	maybe the multiple b_i's are from 1975 Miyamoto, Nagai eqn A1
	where you just evaluate the equation over and over again wich each b_i and sum the functions up?
	but in that case it looks like the multiple b's come with multiple a's and M's?

	--]]
	(function()
		-- where in the [-1,1] unit cube to put the boundary of r=1 of whatever units you use
		-- ... for a sphere ... and this is an ellipsoid
		local coordRadius = .5

		local ngc1560_a = 7.19 * 1e+3 * constants.pc_in_m		-- semi-major axis
		local ngc1560_b = 0.567 * 1e+3 * constants.pc_in_m		-- semi-minor axis
		local ngc1560_M = 1.52e+10 * constants.SolarMass_in_kg	-- total mass

		-- coordinate radius of each axis (used for volume calc)
		local coordRadius_a = coordRadius * (ngc1560_a / ngc1560_a)
		local coordRadius_b = coordRadius * (ngc1560_b / ngc1560_a)

		-- radius .5, grid = 2 * ngc1560_a, so the ellipsoid is ngc1560
		local meter = ngc1560_a / coordRadius

		-- ellipsoid volume in meters:

		-- With this, setting rho=1 will cause the total mass of the grid to be the mass of the galaxy.
		-- 4/3 pi a b^2 rho = M <=> rho = 3 M / (4 pi a b^2)
		local ellipsoidCoordVolume = (4/3) * math.pi * coordRadius_a * coordRadius_b * coordRadius_b
		local kilogram = ngc1560_M / ellipsoidCoordVolume

		local second = 1

		local mins, maxs = createMinsMaxs(coordRadius)

		return {
			name = 'self-gravitation - NGC 1560',

			mins = mins,
			maxs = maxs,
			guiVars = {
				{name = 'a', value = ngc1560_a, units = 'm'},
				{name = 'b', value = ngc1560_b, units = 'm'},
				{name = 'M', value = ngc1560_M, units = 'kg'},
			},
			solverVars = {
				meter = meter,
				kilogram = kilogram,
				second = second,

				-- mind you these will only work if you're using:
				-- 	eqn=euler-lingr
				-- otherwise, these won't be solver_t vars,
				-- and you'll get an error that they're not found
				speedOfLight = constants.speedOfLight_in_m_per_s,
				divPsiWavespeed_g = constants.speedOfLight_in_m_per_s,
				divPhiWavespeed_g = constants.speedOfLight_in_m_per_s,
				gravitationalConstant = constants.gravitationalConstant_in_m3_per_kg_s2,
				coulombConstant = constants.CoulombConstant_in_kg_m3_per_C2_s2,
			},
			--[[
			Now assuming A, B, R, Z in eqn C.1 are in units of [m],
			and M is in units of [kg]
			then rho is in units of [kg/m^3], which is what a density should be.
			But when we get to the normalized-rho def in C.4, we get the result in [a]^5/[a]^8
			for whatever units a,b,r,z are in.
			I'm guessing unitless as well, but then what is their relation (in m) with A,B,R,Z?
			After eqn 4.5 it says r = R/R_0 and z = Z/R_0.
			It doesn't look like a or b are ever defined, but I'm assuming it is done the same way?

			Ok further update: This paper is a mess.  the a,b,lambda are all defined in two places, Section 7 and Appendix D, sometimes contradicting.
			The two sections describe two different profiles of NGC 1560, and some graphs of one model require parameters of the other ... mess of a paper.
			Nowhere is the velocity profile function mentioned.
			It's probably in one of the cited papers, I haven't found the right one yet.
			I have found more complete info between 2021 Ludwig and its cited paper 2006 Cooperstock et al for NGC 3198, so maybe I'll move on to that ...
			--]]
			getInitCondCode = function(self)
				return self.solver.eqn:template[[
	// 2021 Ludwig eqn C.1 ... the non-normalized density function
//// MODULE_DEPENDS: <?=coordMap?>
	real3 const xc = coordMap(x);
	real const rSq = xc.x*xc.x + xc.y*xc.y;	//cylindrical 'r'
	real const a = initCond.a / unit_m;
	real const b = initCond.b / unit_m;
	real const bSq = b*b;
	real const M = initCond.M / unit_kg;
	real const z = xc.z;
	real const zSq = z*z;
	real const bzLenSq = bSq + zSq;
	real const bzLen = sqrt(bzLenSq);
	real const a_plus_bzLen = a + bzLen;
	real const a_plus_bzLen_sq = a_plus_bzLen * a_plus_bzLen;
	real const rSq_plus__a_plus_bzLen_sq = rSq + a_plus_bzLen_sq;
	if (bzLenSq >= 0. && rSq_plus__a_plus_bzLen_sq >= 0.) {
		real const tmp = sqrt(rSq_plus__a_plus_bzLen_sq );
		rho = M / (4. * M_PI
//			* R_0 * R_0 * R_0 * rho_0		// don't use normalized units
		) * bSq * (a * rSq + (a + 3. * bzLen) * a_plus_bzLen_sq)
		/ (
			tmp * tmp * tmp * tmp * tmp
			* bzLen * bzLen * bzLen
		);
	} else {
		rho = 1e-3;
	}

	// ram pressure?
	P = rho * v.dot(v);
	//P = max(P, 1e-7);	//but really, is there any harm in P=0 with the Euler equations?
]]
			end,
		}
	end)(),

	(function()
		local coordRadius = .5
		local ngc3198_rmax = 31.7 * 1e+3 * constants.pc_in_m
		local ngc3198_a = 9.10 * 1e+3 * constants.pc_in_m
		local ngc3198_b = 2.64 * 1e+3 * constants.pc_in_m
		local ngc3198_M = 1.25e+11 * constants.SolarMass_in_kg
		local coordRadius_a = coordRadius * (ngc3198_a / ngc3198_a)
		local coordRadius_b = coordRadius * (ngc3198_b / ngc3198_a)
		local meter = ngc3198_rmax / coordRadius
		local ellipsoidCoordVolume = (4/3) * math.pi * coordRadius_a * coordRadius_b * coordRadius_b
		local kilogram = ngc3198_M / ellipsoidCoordVolume
		local second = 1
		local mins, maxs = createMinsMaxs(coordRadius)
		return {
			name = 'self-gravitation - NGC 3198',
			mins = mins,
			maxs = maxs,
			guiVars = {
				{name = 'a', value = ngc3198_a, units = 'm'},
				{name = 'b', value = ngc3198_b, units = 'm'},
				{name = 'M', value = ngc3198_M, units = 'kg'},
			},
			solverVars = {
				meter = meter,
				kilogram = kilogram,
				second = second,

				-- mind you these will only work if you're using:
				-- 	eqn=euler-lingr
				-- otherwise, these won't be solver_t vars,
				-- and you'll get an error that they're not found
				speedOfLight = constants.speedOfLight_in_m_per_s,
				divPsiWavespeed_g = constants.speedOfLight_in_m_per_s,
				divPhiWavespeed_g = constants.speedOfLight_in_m_per_s,
				gravitationalConstant = constants.gravitationalConstant_in_m3_per_kg_s2,
				coulombConstant = constants.CoulombConstant_in_kg_m3_per_C2_s2,
			},
			getInitCondCode = function(self)
				return self.solver.eqn:template([[
//// MODULE_DEPENDS: <?=coordMap?>
	real3 const xc = coordMap(x);
	real const rSq = xc.x*xc.x + xc.y*xc.y;	//cylindrical 'r'
	real const r = sqrt(rSq);
	real const z = xc.z;

	real const r_in_kpc = r * <?=clnumber(meter / (1e+3 * constants.pc_in_m))?>;
	real const z_in_kpc = z * <?=clnumber(meter / (1e+3 * constants.pc_in_m))?>;

//// MODULE_DEPENDS: Bessel
	// velocity from 2006 Cooperstock et al to fit 1989 Begeman
	real const vmag_per_c = 0
<?
	for _,coeff in ipairs{
		{0.00093352334660, 0.07515079869},
		{0.00020761839560, 0.17250244090},
		{0.00022878035710, 0.27042899730},
		{0.00009325578799, 0.3684854512},
		{0.00007945062639, 0.4665911784},
		{0.00006081834319, 0.5647207491},
		{0.00003242780880, 0.6628636447},
		{0.00003006457058, 0.7610147353},
		{0.00001687931928, 0.8591712228},
		{0.00003651365250, 0.9573314522},
	} do
		local neg_Cn_kn, kn = table.unpack(coeff)
?>			+ <?=neg_Cn_kn?> * BESSJ1(<?=kn?> * r_in_kpc)
<?	end
?>	;

	real const vmag = vmag_per_c * solver->speedOfLight / unit_m_per_s;
	v = real3_real_mul(real3(-xc.y, xc.x, 0.), vmag / r);

	// TODO init GEM potential from 2021 Ludwig once you figure out the graph


	// normrho(r,z=0) for 2021 Ludwig eqn 8.4.b

	real const rspiral = 4.0;
	real const kspiral = 0.1;
	real const gamma0 = 0.95;
	real const gammai = gamma0;
	real const y0 = 8.0;

	real const d = 9200.0;		// kpc

	real const alpha0 = 154.0;
	real const reff = 1.00;
	real const r0 = 6.87;
	real const s0 = 0.586;

<?
local b = {
	0.050231295947568,	-- (from 0.0499)
	-0.000433,
	5.86e-08,
	3.29e-09,
	-1.06e-11,
	1.52e-13,
	2.9e-15,
	-1.75e-17,
}
local d = {
	nil,					-- d1 doesn't exist
	1.9517551593474e-05,	-- (from 1.81e-05)
	-4.96e-07,
	1.85e-09,
	1.07e-11,
	2.04e-14,
	-1.75e-16,
	-1.2e-18,
}
?>

	real const alphae = 316.8;
	real const se = 1.49;

	real Y = 1;
	<? for i=0,4 do ?>{
		real const theta = 2 * M_PI * <?=clnumber(i)?> * kspiral;
		real const ri = rspiral * exp(theta);
		real const yi = y0 * exp(theta * vmag);
		real const dr = r - ri;
		real const gammaisq = gammai * gammai;
		Y += (yi * gammai / M_PI) / (dr * dr + gammaisq);
	}<? end ?>

	real alpha = (180 * 3600 / M_PI) * r / d;

	real s = 0.;
	if (alpha <= alpha0) {
<? for i=8,1,-1 do ?>
		s += <?=clnumber(b[i])?>;
		s *= alpha;
<? end ?>
		s += s0;
	} else if (alpha <= alphae) {
		real const dalpha = alphae - alpha;
<? for i=8,2,-1 do ?>
		s += <?=clnumber(d[i])?>;
		s *= dalpha;
<? end ?>
		s *= dalpha + se;
	} else {
		s = se;
	}

	rho = Y * exp(-pow(r / reff, 1. / s));

	// eqn 6.2: rho(r,z) = rho(r,0) * exp(-z^2 / (2 delta(r)^2 )
	// delta(r) is ... a mess to calculate
	// how about approximating it?
<?
-- piecewise quadratic control points
local pts = {
	{0					,	7.1493036342864			},
	{4.8861861861862	,	6.7384253718743			},
	{6.852952952953		,	6.361204928482			},
	{8.8197197197197	,	5.8907000621838			},
	{9.8031031031031	,	5.6296706252255			},
	{10.786486486486	,	5.3571550559966			},
	{11.76986986987		,	5.0775983990215			},
	{12.999099099099	,	4.7247306823609			},
	{13.982482482482	,	4.444299392457			},
	{14.965865865866	,	4.1688493246293			},
	{15.949249249249	,	3.9006182904375			},
	{17.916016016016	,	3.3908327257272			},
	{19.882782782783	,	2.918292334253			},
	{21.84954954955		,	2.4781486437422			},
	{23.816316316316	,	2.0600568093525			},
	{25.783083083083	,	1.6481068418665			},
	{27.74984984985		,	1.214040418122			},
	{28.733233233233	,	0.96898392138409		},
	{29.716616616617	,	0.67271150250128		},
	{30.454154154154	,	0.34095094607551		},
	{30.7				,	0.095697556491137		},
}
for i=1,#pts-1,2 do
?>	<? if i > 1 then ?>} else <? end ?>if (r_in_kpc < <?=clnumber(pts[i+2][1])?>) {
<?
	local x1, y1 = table.unpack(pts[i])
	local x2, y2 = table.unpack(pts[i+1])
	local x3, y3 = table.unpack(pts[i+2])
	local x0 = x1
	x3 = x3 - x0
	x2 = x2 - x0
	x1 = x1 - x0
	-- https://math.stackexchange.com/a/680695/206369
	local a = (x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1)) / ((x1 - x2) * (x1 - x3) * (x2 - x3))
	local b = (y2 - y1) / (x2 - x1) - a * (x1 + x2)
	local c = y1 - a * x1 * x1 - b * x1
?>		real const delta_for_r = <?=clnumber(c)?> + (r_in_kpc - <?=clnumber(x0)?>) * (<?=clnumber(b)?> + (r_in_kpc - <?=clnumber(x0)?>) * <?=clnumber(a)?>);
		real const z_delta_ratio = z_in_kpc / delta_for_r;
		real const zinfl = exp(-.5 * z_delta_ratio * z_delta_ratio);
		rho *= zinfl;
<?
end
?>
	} else {
		// out of bounds density
		rho = 0.;
	}

	//ok now assume rho is in [0,1] ... so scale velocity by rho/rhomax
	// but in reality (also in the paper) the graph of normrho is only [0,1] for the data sample points
	// for the Ludwig function fitting it is about [0, 1.3]
	real const normrhomax = 1.3587387933273;	//sup of numerics graph of normrho(r)
	v = real3_real_mul(v, rho / normrhomax);		// assuming rho is normalized at this point

	//then inc rho past some min value (since I don't handle vacuum cells atm)
	real const rhomin = 1e-4;	// this is numeric / used to avoid vacuum
	rho += rhomin;

	// at this point I'm assuming the normalized units are close enough
	// but we can be certain:
	real const rho0 = 6.54e-21;	// in kg/m^3 ... this is from Ludwig 2021 after eqn 8.5
	rho *= rho0 / unit_kg_per_m3;
	// hmm, after this rho is about [0,10] in the simulation units ... which means my estimation of unit_kg above based on the galactic properties is about 10x off ...


	// galactic pressure ... ram pressure?
	P = rho * v.dot(v);
	//P = max(P, 1e-7);	//but really, is there any harm in P=0 with the Euler equations?
	//P *= zinfl;
]],				{
					constants = constants,
					meter = meter,
				})
			end,
		}
	end)(),

	SelfGravProblem:subclass{
		name = 'self-gravitation test 1',
		-- hmm, what about spherical coordinates...
		--mins = {-1,-1,-1},
		--maxs = {1,1,1},
		sources = {
			{center={0, 0, 0}, radius = .5},
		},
		getRadiusCode = function(source)
			return '.5 - .01 * (U.s[0] - .5)'
		end,
	},

	SelfGravProblem:subclass{
		name = 'self-gravitation test 1 spinning',
		sources={
			{
				center={0, 0, 0},
				radius = .2,
				-- srhd solver requires the max velocity not to exceed 1 ...
				inside = [[
	v.x = -2 * delta.y;
	v.y = 2 * delta.x;
	rho = 1.;
	P = 1.;
]],
			},
		},
		getRadiusCode = function(source)
			-- TODO compute dr's component in each dx^i's, and scale our random number by that
			-- U.rho holds noise
			return '.2 - .005 * (U.s[0] - .5)'
		end,
	},

	SelfGravProblem:subclass{
		name = 'self-gravitation test 2',
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

	SelfGravProblem:subclass{
		-- TODO add tidal-locked rotations
		name = 'self-gravitation test 2 orbiting',
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

	SelfGravProblem:subclass{
		name = 'self-gravitation test 4',
		sources={
			{center={.25, .25, 0}, radius = .1},
			{center={-.25, .25, 0}, radius = .1},
			{center={.25, -.25, 0}, radius = .1},
			{center={-.25, -.25, 0}, radius = .1},
		},
	},

	{
		name = 'self-gravitation soup',
		getInitCondCode = function(self)
			return [[
	rho = .1 * U.s[0] + .1;
	v.x = .2 * U.s[1] - .1;
	v.y = .2 * U.s[2] - .1;
	v.z = .2 * U.s[3] - .1;
	P = .1 * U.s[4] + .1;
]]
		end,
	},

	SelfGravProblem:subclass{
		name = 'self-gravitation Jeans, right?',
		getRadiusCode = function(source)
			return '.5 - .02 * (U.s[0] - .5)'
		end,
		sources = {
			{center={0, 0, 0}, radius = .5, inside=[[
	P = 1;
	rho = .1;
	if (distSq > .4 * .4) {
		rho = 1.;
	}
]]},
		},
	},



	{
		name = 'Maxwell default',
		solverVars = {
			--meter = constants.speedOfLight_in_m_per_s,
			--speedOfLight = constants.speedOfLight_in_m_per_s,
		},
		getInitCondCode = function(self)
			return self.solver.eqn:template([[
	D.x = <?=scalar?>(lhs ? 1 : -1);
	D.y = <?=scalar?>(1);
	B.y = <?=scalar?>(-1);
	B.z = <?=scalar?>(lhs ? 1 : -1);
]])
		end,
	},

	-- now that I think about it, this is dumb, because the derivative will be zero.
	-- oscillating boundary is much better.
	{
		name = 'Maxwell constant',
		guiVars = {
			{name = 'Dx', value = 1},
			{name = 'Dy', value = 0},
			{name = 'Dz', value = 0},
			{name = 'Bx', value = 1},
			{name = 'By', value = 0},
			{name = 'Bz', value = 0},
		},
		getInitCondCode = function(self)
			return self.solver.eqn:template([[
	D.x = <?=scalar?>(initCond.Dx);
	D.y = <?=scalar?>(initCond.Dy);
	D.z = <?=scalar?>(initCond.Dz);
	B.x = <?=scalar?>(initCond.Bx);
	B.y = <?=scalar?>(initCond.By);
	B.z = <?=scalar?>(initCond.Bz);
]])
		end,
	},

	{
		name = 'Maxwell empty waves',
		getInitCondCode = function(self)
			local solver = assert(self.solver)
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
		name = 'Maxwell empty waves',
		getInitCondCode = function(self)
			local solver = assert(self.solver)
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
		getInitCondCode = function(self)
			addMaxwellOscillatingBoundary{
				solver = self.solver,
				side = 'xmin',
				dir = 'x',
				amplitude = 1,
				period = 10,
			}
			return self.solver.eqn:template[[
//// MODULE_DEPENDS: <?=coordMap?>
	real3 const xc = coordMap(x);
	if (real3_lenSq(xc) < .2*.2) {
		//2018 Balezin et al "Electromagnetic properties of the Great Pyramids..."
//// MODULE_DEPENDS: cplx
		permittivity = <?=susc_t?>(cplx(5., .1));
	}
]]
		end,
	},

	{
		name = 'Maxwell scattering around pyramid',
		getCodePrefix = function(self)
			return self.solver.eqn:template([[
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
	local p = table(pn.p):mapi(clnumber):concat', '
	local n = table(pn.n):mapi(clnumber):concat', '
?>
			&& (xc - real3(<?=p?>).dot(real3(<?=n?>)) < 0.
<? end ?>
		)
<? end ?>
	;
}
]])
		end,
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			-- hmm, choosing min or max doesn't matter, it always shows up on min...
			addMaxwellOscillatingBoundary{
				solver = solver,
				side = xNames[solver.dim]..'max',
				dir = xNames[solver.dim],
				amplitude = -1,
				period = 10,
			}
			return solver.eqn:template([[
//// MODULE_DEPENDS: <?=coordMap?>
	real3 xc = coordMap(x);
	xc = real3_real_mul(xc, 2.);
	if (testTriangle(xc)) {
		//2018 Balezin et al "Electromagnetic properties of the Great Pyramids..."
//// MODULE_DEPENDS: cplx
		permittivity = <?=susc_t?>(cplx(5., .1));
	}
]])
		end,
	},

	{
		name = 'Maxwell scattering around square',
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			-- hmm, choosing min or max doesn't matter, it always shows up on min...
			addMaxwellOscillatingBoundary{
				solver = solver,
				side = xNames[solver.dim]..'max',
				dir = xNames[solver.dim],
				amplitude = -1,
				period = 10,
			}
			return solver.eqn:template([[
//// MODULE_DEPENDS: <?=coordMap?>
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
//// MODULE_DEPENDS: cplx
		permittivity = <?=susc_t?>(cplx(5., .1));
	}
]])
		end,
	},

	{
		name = 'Maxwell scattering around Koch snowflake',
		getCodePrefix = function(self)
			return self.solver.eqn:template([[

#define sqrt3 <?=clnumber(math.sqrt(3))?>

#define p1 real3(0, .5, 0)
#define p2 real3(.25*sqrt3, -.25, 0.)
#define p3 real3(-.25*sqrt3, -.25, 0.)

#define n1 real3(0, 1, 0)
#define n2 real3(.5*sqrt3, -.5, 0)
#define n3 real3(-.5*sqrt3, -.5, 0)

//initial branches
<? for i=1,3 do ?>
real3 branch<?=i?>(real3 x) {
	x = real3(-x.x, -x.y, 0);	//180 rotation
	real3 n = real3(-n<?=i?>.y, n<?=i?>.x, 0);	//angle of rotation of the normal
	x = real3(n.x * x.x - n.y * x.y, n.x * x.y + n.y * x.x, 0.);	//rotate by 'n'
	x.y -= sqrt3*3./8.;	//translate to center
	x = real3_real_mul(x, 3.);	//scale up by 3
	x = real3(-x.x, -x.y, 0);	//180 rotation
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
	x = real3(x.x*c - x.y*s, x.x*s + x.y*c, 0.);	//rotate by c,s
	x = real3_real_mul(x, 3);	//scale by 3
	x.y += sqrt3;	//translate to center
	return x;
}

real3 branch2_4(real3 x) {
	real c = .5;
	real s = -sqrt3*.5;
	x = real3(x.x*c - x.y*s, x.x*s + x.y*c, 0.);	//rotate by c,s
	x = real3_real_mul(x, 3);	//scale by 3
	x.y += sqrt3;	//translate to center
	return x;
}

bool testTriangle(real3 xc) {
	return ((xc - p1).dot(n1) < 0. &&
		(xc - p2).dot(n2) < 0. &&
		(xc - p3).dot(n3) < 0.);
}
]])
		end,
		getInitCondCode = function(self)
			local solver = assert(self.solver)
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

			return solver.eqn:template([[
//// MODULE_DEPENDS: <?=coordMap?>
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
//// MODULE_DEPENDS: cplx
		permittivity = <?=susc_t?>(cplx(5., .1));
	}

]], 		{
				resistivities = resistivities,
			})
		end,
	},

	{
		name = 'Maxwell Lichtenberg',
		getInitCondCode = function(self)
			local solver = assert(self.solver)
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
					solver.modules:getCodeAndHeader(solver.sharedModulesEnabled:keys():unpack())
						:gsub('//// BEGIN INCLUDE FOR FFI_CDEF.-//// END INCLUDE FOR FFI_CDEF', '')
					,
					solver.eqn:template([[
//single cell domain
kernel void addExtraSource(
	global <?=cons_t?>* UBuf
) {
	UBuf[INDEX(solver, <?=src[1]?>,<?=src[2]?>,<?=src[3]?>)].D.x = -10;
	UBuf[INDEX(solver, <?=dst[1]?>,<?=dst[2]?>,<?=dst[3]?>)].D.x = -10;
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
		resetState = function(self)
			local solver = assert(self.solver)
			-- super calls applyInitCondKernel ...
			EulerInitCond.resetState(self)
			-- and here I'm going to fill the permittivity 'eps' with random noise
			-- ... and put a source + and - current 'sigma' at two points on the image
			local ptr = ffi.cast(solver.eqn.symbols.cons_t..'*', solver.UBufObj:toCPU())
			for i=0,solver.numCells-1 do
				ptr[i].sigma = math.random() * 1e-4 + 1e-7
			end
			solver.UBufObj:fromCPU(ffi.cast('real*', ptr))
		end,
	},

	{
		name = 'Maxwell wire',
		getInitCondCode = function(self)
			local solver = assert(self.solver)
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
			}:mapi(function(v) return v * Ohm_in_m end)
			return solver.eqn:template([[
	D.x = <?=scalar?>(1.);

	//conductivity = <?=susc_t?>(<?=clnumber(1/resistivities.air)?>);

	real r2 = x.y * x.y<? if solver.dim == 3 then ?> + x.z * x.z<? end ?>;

	if (r2 < .1*.1) {
		//conductivity = <?=susc_t?>(<?=clnumber(1/resistivities.copper)?>);
		permittivity = <?=susc_t?>(cplx(5., .1));
#error TODO assign rhoCharge and J
	}
]], 		{
				resistivities = resistivities,
			})
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
		getInitCondCode = function(self)
			return self.solver.eqn:template[[
	D.z = <?=scalar?>( E0 * sin(m * M_PI * x.x / x0) * sin(n * M_PI * x.y / y0) );
]]
		end,
	},

	{
		guiVars = {
			{name = 'rhoCharge0', value = 1},
		},
		name = 'Maxwell charged particle',
		getInitCondCode = function(self)
			return self.solver.eqn:template[[
	int4 const i = globalInt4();
	rhoCharge = (i.x == solver->gridSize.x/2 && i.y == solver->gridSize.y/2 && i.z == solver->gridSize.z/2) ? initCond.rhoCharge0 : 0.;
]]
		end,
	},

	{
		name = '2017 Degris et al',
		getInitCondCode = function(self)
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
		getInitCondCode = function(self)
			-- I am not correctly modeling the top boundary
			self.solver:setBoundaryMethods{
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
		getInitCondCode = function(self)
			return self.solver.eqn:template[[
//// MODULE_DEPENDS: <?=coordMap?>
	real3 xc = coordMap(x);
	real r = xc.length();
	P = initCond.P;
	rho = initCond.rho;
	v.x = -xc.y * initCond.v / r;
	v.y = xc.x * initCond.v / r;
	real s = sign(r - .5);
	B.x = -xc.y * s * initCond.B / r;
	B.y = xc.x * s * initCond.B / r;
]]
		end,
	},

	{
		name = 'two-fluid EMHD soliton ion',
		getInitCondCode = function(self)
			return [[
	real const L = 12.;
	rho = 1. + (real)exp((real)-25. * (real)fabs(x.x - L / (real)3.));
]]
		end,
	},
	{
		name = 'two-fluid EMHD soliton electron',
		getInitCondCode = function(self)
			return [[
	real const L = 12.;
	rho = 5. * (1. + (real)exp((real)-25. * (real)fabs(x.x - L / (real)3.)));
]]
		end,
	},
	{
		name = 'two-fluid EMHD soliton maxwell',
		getInitCondCode = function(self)
			error("TODO")
		end,
	},

	--1999 Hudson - "Numerical Techniques for the Shallow Water Equations"
	-- TODO ... this looks like problem B ... but I have problem B below ... hmm ...
	-- the second half looks dif tho ... where is it from?
	{
		name = 'shallow water constant',
		solverVars = {
			water_D = .5,
		},
		getInitCondCode = function(self)
			return [[
	real const s = .5 * (x.x + 1.);

	//this is bathymetry depth
	real water_B = 0.;
#if 0
	if (.4 <= s && s <= .6) {
		water_B = (1./8.) * (cos(10. * M_PI * (s - .5)) + 1.);
	}
#else
/*
	with g = 1, D = 0, with source disabled:

	hll: mx,t : 0.05 -> 0.045	<- this is the correct flux value
	when roeUseFluxFromCons==false:
	roe: mx,t : 0.1 -> 0.09		<- this is doubled, due to dF/dU * U doubling the g h^2 term compared to F's 1/2 g h^2
	when roeUseFluxFromCons==true:	<- this is the correct flux value
	hll: mx,t : 0.05 -> 0.045

	hll: Fh : .0000105 -> .00000905
	roe: Fh : .0000105 -> .00000905

	hll: Fmx : 0.5 -> 0.405
	when roeUseFluxFromCons==false:
	roe: Fmx : 1 -> 0.81

... why is the roe derivative half of the hll derivative?
is it always , in all cases?

what should this be?
g = 1
U = [h, hv]
F = [hv, hv^2 + .5 g h^2]
dF/dx = [
	(v) dh/dx + (h) dv/dx,
	(v^2 + g h) dh/dx + (h v) dv/dx
]
h = 1 - .1 s, s in [0,1], s = .5 (x + 1)
v = 0
U = [1 - .1 s, 0]
F = [0, .5 (1 - .1 s)^2]
	= [0, .5  to .405]
dF/dx = [0, .1 g (1 - .1 s)]
	  = [0, .1  - .01 s]
	  = from .1 to .09
So roe is giving dF/dx, hll is giving F
Why do I get the feeling that in some places I'm using F in place of dF/dx ?
Turns out this mixup is used throughout literature
because it just so happens that, for the Euler fluid equations, dF/dU * U = F, so everyone assumes it's true.
Probably because the U terms within dF/dU happen to be linear?
Whereas for the shallow-water equation one term is quadratic (which makes the 1/2 disappear in the derivative).

So new question, which is correct?
Do we want 'calcFluxForInterface' to calculate F, which fvsolver then linearly approximates as dF/dx ~ (FR-FL)/(xR-xL)?
Or likewise by the divergence theorem that int div F dV = int F dot n dS?
In both cases it looks like F is wanted, not dF/dU.
*/
	water_B = .1 * s;
#endif

	real const water_H = solver->water_D - water_B;
	depth = water_H;

	// this is wave height for shallow water equations:
	real const water_h = 1. - water_B;

	//here's our placeholder variable I call 'rho' just for compat with euler fluid equation code
	rho = water_h;
]]
		end,
	},

	--1999 Hudson - "Numerical Techniques for the Shallow Water Equations"
	-- problem A
	{
		name = 'shallow water problem A',
		mins = {0,0,0},
		maxs = {1,1,1},
		solverVars = {
			water_D = .5,
		},
		getInitCondCode = function(self)
			return [[
	real const s = (x.x - solver->initCondMins.x) / (solver->initCondMaxs.x - solver->initCondMins.x);

	real water_B = 0;	//page4: "This is due to the riverbed being of constant depth..."

	real const phi0 = 0.5;	//only really specified on page 5, not said to be all phi0's, not said why ...
	real water_h;
	if (s < .5) {
		water_h = 1. - water_B;
	} else {
		water_h = phi0 - water_B;
	}
	rho = water_h;
]]
		end,
	},

	--1999 Hudson - "Numerical Techniques for the Shallow Water Equations"
	-- problem B
	{
		name = 'shallow water problem B',
		mins = {0,0,0},
		maxs = {1,1,1},
		guiVars = {
			{name = 'phi0', value = .5},
		},
		solverVars = {
			water_D = .5,
		},
		getInitCondCode = function(self)
			return [[
	real const s = (x.x - solver->initCondMins.x) / (solver->initCondMaxs.x - solver->initCondMins.x);

	//this is bathymetry depth
	real water_B = 0.;
	if (.4 <= s && s <= .6) {
		water_B = (1./8.) * (cos(10. * M_PI * (s - .5)) + 1.);
	}

	real const water_H = solver->water_D - water_B;
	depth = water_H;

	// this is wave height for shallow water equations:
	real water_h;
	if (s < .5) {
		water_h = 1 - water_B;
	} else {
		water_h = initCond.phi0 - water_B;
	}
	//here's our placeholder variable I call 'rho' just for compat with euler fluid equation code
	rho = water_h;
]]
		end,
	},

	--1999 Hudson - "Numerical Techniques for the Shallow Water Equations"
	-- problem C
	{
		name = 'shallow water problem C',
		solverVars = {
			water_D = .5,
		},
		getInitCondCode = function(self)
			return [[
	real const s = (x.x - solver->initCondMins.x) / (solver->initCondMaxs.x - solver->initCondMins.x);

	//this is bathymetry depth
	real water_B = 0.;
	if (.4 <= s && s <= .6) {
		water_B = (1./8.) * (cos(10. * M_PI * (s - .5)) + 1.);
	}

	real const water_H = solver->water_D - water_B;
	depth = water_H;

	// this is wave height for shallow water equations:
	real water_h;
	if (s < .1) {
		water_h = 1. - water_B;
	} else if (s < .2) {
		water_h = 1.2 - water_B;
	} else {
		water_h = 1. - water_B;
	}

	//here's our placeholder variable I call 'rho' just for compat with euler fluid equation code
	rho = water_h;
]]
		end,
	},

	{
		name = 'shallow water parabola',
		getInitCondCode = function(self)
			return [[
	real const s = (x.x - solver->initCondMins.x) / (solver->initCondMaxs.x - solver->initCondMins.x);

	//this is bathymetry depth
	real const water_B = 2. * x.x * x.x - 1.;

	real const water_H = solver->water_D - water_B;
	depth = water_H;

	// this is wave height for shallow water equations:
	real water_h;
	if (s < .1) {
		water_h = -water_B;
	} else if (s < .2) {
		water_h = .2 - water_B;
	} else {
		water_h = -water_B;
	}

	if (water_h < water_H) water_h = water_H;
	if (water_h < 0) water_h = 0;

	//here's our placeholder variable I call 'rho' just for compat with euler fluid equation code
	rho = water_h;
]]
		end,
	},


	-- 2003 Rogers et al "Mathematical balancing of flux gradient and source terms prior to using Roe's approximate Riemann solver"
	-- section 4.2.3
	-- TODO boundary condition with inflow of 0.18 m^2 / s (after renomralizing my domain from [-1,1] to [0,25])
	-- and then depth at right boundary is 0.33 m
	-- ... what's the depth at left boundary?
	-- Figure 4.c looks like it's constant acceleration across the grid of 0.18 m/s^2 (m^2/s a typo?)
	-- still, what's the left boundary depth?
	{
		name = '2003 Rogers',
		getInitCondCode = function(self)
			return [[
	real const s = (x.x - solver->initCondMins.x) / (solver->initCondMaxs.x - solver->initCondMins.x);

	if (8./25. < s && s < 12./25.) {
		water_B = 0.05 * (s - 10./25.) * (s - 10./25.);
	} else {
		water_B = 0.2;
	}

	real const water_H = solver->water_D - water_B;
	depth = water_H;

	rho = -water_B;
]]
		end,
	},

}:mapi(function(cl)
	if EulerInitCond:isa(cl) then
		return cl:subclass()
	end
	return EulerInitCond:subclass(cl)
end)

function EulerInitCond:getList()
	return initConds
end

return EulerInitCond

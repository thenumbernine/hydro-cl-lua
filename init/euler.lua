local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local clnumber = require 'cl.obj.number'
local template = require 'template'
local ffi = require 'ffi'
local InitState = require 'init.init'

local function quadrantProblem(args)
	args.initState = function(self, solver)
		solver.cfl = .475
		solver:setBoundaryMethods'freeflow'
		local function build(i)
			local q = args[i]
			return table.map(q, function(v,k,t)
				return k..'='..v..';', #t+1
			end):concat' '
		end
		return [[
	bool xp = x.x > mids.x;
	bool yp = x.y > mids.y;
	if (yp) {
		if (xp) {
			]]..build(1)..[[
		} else {
			]]..build(2)..[[	
		}
	} else {
		if (!xp) {
			]]..build(3)..[[
		} else {
			]]..build(4)..[[
		}
	}
]]
	end
	return args
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

function SelfGravProblem:__call(solver)
	local args = self.args
	solver.useGravity = true
	--[[ the boundary might not be necessary/appropriate, esp for cylindrical geometry
	solver:setBoundaryMethods'freeflow'
	--]]

	return template([[
	rho = .1;
	P = 1;
	//v[i] = .1 * noise * crand();
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

local function addMaxwellOscillatingBoundary(solver)
	-- this args is only for the UBuf boundary program -- not calle for the Poisson boundary program
	function solver:getBoundaryProgramArgs()
		-- i'm completely overriding this
		-- so I can completely override boundaryMethods for the solver boundary kernel
		-- yet not for the poisson boundary kernel
		local boundaryMethods = table(self.boundaryMethods)
		boundaryMethods.xmin = function(args)
			local U = 'buf['..args.index'j'..']'
			return template([[
	<?=U?>.epsE.y = (real)sin((real)10. * t) / <?=U?>.eps;
]], {U=U})
		end

		-- same as super 
		-- except with extraAgs
		-- and using boundaryMethods instead of self.boundaryMethods
		-- (which I leave alone so Poisson can deal with it)
		return {
			type = self.eqn.cons_t,
			extraArgs = {'real t'},
			-- remap from enum/combobox int values to functions from the solver.boundaryOptions table
			methods = table.map(boundaryMethods, function(v)
				if type(v) == 'function' then return v end
				return (select(2, next(self.boundaryOptions[v])))
			end),
			mirrorVars = self.eqn.mirrorVars,
		}
	end
	
	-- this runs before refreshBoundaryProgram, so lets hijack refreshBoundaryProgram and add in our time-based boundary conditions
	local oldBoundary = solver.boundary
	function solver:boundary()
		self.boundaryKernel:setArg(1, ffi.new('real[1]', self.t))
		oldBoundary(self)
	end
end

local initStates = table{
	{
		name = 'constant',
		initState = function(self, solver)
			if solver.eqn.guiVars.heatCapacityRatio then	
				solver.eqn.guiVars.heatCapacityRatio.value = 7/5
			end
			return [[
	rho = 1;
	P = 1;
]]
		end,
	},
	{
		name = 'linear',
		initState = function(self, solver)
			if solver.eqn.guiVars.heatCapacityRatio then	
				solver.eqn.guiVars.heatCapacityRatio.value = 7/5
			end
			return [[
	rho = 1 + x.x;
	P = 1 + x.x;
]]
		end,
	},
	{
		name = 'constant with motion',
		initState = function(self, solver) 
			return [[
	rho=1;
	v.x=1;
	v.y=1;
	v.z=1;
	P=1;
]]
		end,
	},
	{
		name = 'linear',
		initState = function(self, solver)
			return '	rho=2+x.x; P=1;'
		end,
	},
	{
		name = 'gaussian',
		initState = function(self, solver)
			return [[
	real sigma = 1. / sqrt(10.);
	real xSq = real3_dot(x,x);
	rho = exp(-xSq / (sigma*sigma)) + .1;
	P = 1 + .1 * (exp(-xSq / (sigma*sigma)) + 1) / ((heatCapacityRatio - 1.) * rho);
]]
		end,
	},
	{
		name = 'advect wave',
		initState = function(self, solver)
			return [[
	real3 xc = real3_sub(coordMap(x), _real3(-.5, 0, 0));
	real rSq = real3_lenSq(xc);
	rho = exp(-100*rSq) + 1.;
	v.x = 1;
	P = 1;
]]
		end,
	},

	
	{
		name = 'Sod',
		initState = function(self, solver)
			if solver.eqn.guiVars.heatCapacityRatio then	
				solver.eqn.guiVars.heatCapacityRatio.value = 7/5
			end
			return [[
	rho = lhs ? 1. : .125;
	P = lhs ? 1. : .1;
]]
		end,
	},
	{
		name = 'Sedov',
		initState = function(self, solver)
			return [[
	rho = 1;
	P = (i.x == gridSize.x/2 && i.y == gridSize.y/2 && i.z == gridSize.z/2) ? 1e+3 : 1;
]]
		end,
	},
	-- http://www.astro.uni-bonn.de/~jmackey/jmac/node7.html
	-- http://www.astro.princeton.edu/~jstone/Athena/tests/brio-wu/Brio-Wu.html
	{
		name = 'Brio-Wu',
		initState = function(self, solver)
			if solver.eqn.guiVars.heatCapacityRatio then	
				solver.eqn.guiVars.heatCapacityRatio.value = 2
			end
			return [[
	rho = lhs ? 1 : .125;
	P = lhs ? 1 : .1;
	B.x = .75;
	B.y = lhs ? 1 : -1;
	B.z = 0;
]]
		end,
	},
	-- http://www.astro.virginia.edu/VITA/ATHENA/ot.html
	-- http://www.astro.princeton.edu/~jstone/Athena/tests/orszag-tang/pagesource.html
	{
		name = 'Orszag-Tang',
		initState = function(self, solver)
			if solver.eqn.guiVars.heatCapacityRatio then	
				solver.eqn.guiVars.heatCapacityRatio.value = 5/3
			end
			return [[
	const real B0 = 1./sqrt(4. * M_PI);
	rho = 25./(36.*M_PI);
	v.x = -sin(2.*M_PI*(x.y+.5));
	v.y = sin(2.*M_PI*(x.x+.5));
	v.z = 0;
	P = 5./(12.*M_PI);	// is this hydro pressure or total pressure?
	B.x = -B0 * sin(2. * M_PI * (x.y+.5));
	B.y = B0 * sin(4. * M_PI * (x.x+.5));
	B.z = 0;
]]
		end,
	},
	-- http://www.cfd-online.com/Wiki/Explosion_test_in_2-D
	{
		name = 'sphere',
		initState = function(self, solver)
			return [[
	real rSq = real3_dot(x,x);
	bool inside = rSq < .5*.5;
	rho = inside ? 1 : .1;
	P = inside ? 1 : .1;
]]
		end,
	},
	{
		name = 'rarefaction wave',
		initState = function(self, solver)
			return [[
	real delta = .1;
	rho = 1;	// lhs ? .2 : .8;
	v.x = lhs ? .5 - delta : .5 + delta;
	P = 1;
]]
		end,
	},
	
	-- 2D tests described in Alexander Kurganov, Eitan Tadmor, Solution of Two-Dimensional Riemann Problems for Gas Dynamics without Riemann Problem Solvers
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
		initState = function(self, solver)
			solver.cfl = .5	-- needs a slower cfl
			if solver.eqn.guiVars.heatCapacityRatio then	
				solver.eqn.guiVars.heatCapacityRatio.value = 4/3
			end
			return [[
	rho = 1;
	v.x = 1. - 1e-5;
	P = (heatCapacityRatio - 1.) * rho * (1e-7 / sqrt(1. - v.x * v.x));
]]
		end,
	},
	{
		name = 'relativistic blast wave test problem 1',
		initState = function(self, solver)
			solver.cfl = .5	-- needs a slower cfl
			if solver.eqn.guiVars.heatCapacityRatio then	
				solver.eqn.guiVars.heatCapacityRatio.value = 5/3
			end
			return [[
	rho = lhs ? 10 : 1;
	P = (heatCapacityRatio - 1.) * rho * (lhs ? 2 : 1e-6);
]]
		end,
	},
	{
		name = 'relativistic blast wave test problem 2',
		initState = function(self, solver)
			solver.cfl = .5	-- needs a slower cfl
			if solver.eqn.guiVars.heatCapacityRatio then	
				solver.eqn.guiVars.heatCapacityRatio.value = 5/3
			end
			return [[
	rho = 1;
	P = lhs ? 1000 : .01;
]]
		end,
	},
	{
		name = 'relativistic blast wave interaction',
		initState = function(self, solver)
			solver.cfl = .5	-- needs a slower cfl
			return template([[
	
	bool wave1 = true
<? for i=0,solver.dim-1 do ?>
		&& x.s<?=i?> < .9 * mins.s<?=i?> + .1 * maxs.s<?=i?>
<? end ?>;
	bool wave2 = true
<? for i=0,solver.dim-1 do ?>
		&& x.s<?=i?> > .1 * mins.s<?=i?> + .9 * maxs.s<?=i?>
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
		initState = function(self, solver)
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
	{
		name = 'Kelvin-Hemholtz',
		initState = function(self, solver)
			solver:setBoundaryMethods'periodic'
			
			local moveAxis = 0
			local sliceAxis = 1
			
			-- move around the cylinder
			if solver.geometry.name == 'cylinder' then
				moveAxis = 1
				sliceAxis = 0
			end
			
			return template([[
	real yq1 = mins.s<?=sliceAxis?> * .75 + maxs.s<?=sliceAxis?> * .25;
	real yq2 = mins.s<?=sliceAxis?> * .25 + maxs.s<?=sliceAxis?> * .75;
	bool inside = (x.s<?=sliceAxis?> > yq1 && x.s<?=sliceAxis?> < yq2);

	real theta = 2. * M_PI;
<?
for i=0,solver.dim-1 do 
	if i ~= sliceAxis then
?>	theta *= (x.s<?=i?> - mins.s<?=i?>) / (maxs.s<?=i?> - mins.s<?=i?>);
<? 
	end	
end ?>

	real noise = (maxs.x - mins.x) * 1e-4;
	rho = inside ? 2 : 1;
	v.x = cos(theta) * noise;
#if dim == 2
	v.y = sin(theta) * noise;
#endif
#if dim == 3
	v.z = sin(theta) * noise;
#endif
	v.s<?=moveAxis?> += (inside ? -.5 : .5);
	P = 2.5;
]],				{
					solver = solver,
					sliceAxis = sliceAxis,
					moveAxis = moveAxis,
				}
			)
		end,
	},
	
	-- http://www.astro.princeton.edu/~jstone/Athena/tests/rt/rt.html
	-- TODO fixme 
	{
		name = 'Rayleigh-Taylor',
		initState = function(self, solver)
			local xs = {'x', 'y', 'z'}
			
			local xmid = (solver.mins + solver.maxs) * .5
			
			-- triple the length along the interface dimension
			local k = solver.dim
			solver.mins[k] = xmid[k] + (solver.mins[k] - xmid[k]) * 3
			solver.maxs[k] = xmid[k] + (solver.maxs[k] - xmid[k]) * 3

			-- triple resolution along interface dimension
			--size[k] = size[k] * 3
			-- (can it handle npo2 sizes?)

			local boundaryMethods = {}
			for i,x in ipairs(xs) do
				for _,minmax in ipairs{'min', 'max'} do
					boundaryMethods[x..minmax] = i == solver.dim and 'mirror' or 'periodic'
				end
			end
			solver:setBoundaryMethods(boundaryMethods)
		
			return [[
	const real3 externalForce = _real3(0,1,0);
	ePot = 0. <? 
for side=0,solver.dim-1 do
?> + (x.s<?=side?> - mins.s<?=side?>) * externalForce.s<?=side?><?
end ?>;
	int topdim = <?=solver.dim-1?>;
	bool top = x.s[topdim] > mids.s[topdim];
	//TODO noise too?
	rho = top ? 2 : 1;
	P = 2.5 - rho * ePot;
	// or maybe it is ... pressure = (gamma - 1) * density * (2.5 - potentialEnergy)
]]
		end,
	},


	--http://www.astro.virginia.edu/VITA/ATHENA/dmr.html
	{
		name = 'double mach reflection',
		initState = function(self, solver)
			-- I am not correctly modeling the top boundary
			solver.mins = vec3(0,0,0)
			solver.maxs = vec3(4,1,1)
			solver:setBoundaryMethods{
				xmin = 'freeflow',
				xmax = 'freeflow',
				ymin = 'mirror',
				ymax = 'freeflow',
				zmin = 'mirror',
				zmax = 'mirror',
			}
			if solver.eqn.guiVars.heatCapacityRatio then	
				solver.eqn.guiVars.heatCapacityRatio.value = 7/5
			end
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
		initState = function(self, solver)
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
		name='shock bubble interaction',
		initState = function(self, solver)
			solver:setBoundaryMethods'freeflow'
			return [[
	const real waveX = -.45;
	real3 bubbleCenter = _real3(0,0,0);
	real bubbleRadius = .2;
	real3 delta = real3_sub(x, bubbleCenter);
	real bubbleRSq = real3_lenSq(delta);
	rho = x.x < waveX ? 1. : (bubbleRSq < bubbleRadius*bubbleRadius ? .1 : 1);
	P = x.x < waveX ? 1 : .1;
	v.x = x.x < waveX ? 0 : -.5;
]]
		end,
	},

	-- gravity potential test - equilibrium - Rayleigh-Taylor ... or is it Jeans? (still has an shock wave ... need to fix initial conditions?)

	{
		name = 'self-gravitation test 1',
		initState = function(self, solver)
			local radius = .5
			local f = SelfGravProblem{
				solver = solver,
				sources={
					{center={0, 0, 0}, radius = radius},
				},
			}
			if solver.geometry.name == 'cylinder' then
				solver:setBoundaryMethods{
					xmin = 'freeflow',
					xmax = 'freeflow',
					ymin = 'periodic',
					ymax = 'periodic',
					zmin = 'freeflow',
					zmax = 'freeflow',
				}
				return template([[
	P = 1;
	rho = x.s[0] < <?=radius?> ? 1 : .1;
]], {radius=radius})
			else
				return f(solver)
			end
		end
	},

	{
		name = 'self-gravitation test 1 spinning',
		initState = function(self, solver)
			local inside = [[
	v.x = -2 * delta.y;
	v.y = 2 * delta.x;
	rho = 1.;
	P = 1.;
]]

			return SelfGravProblem{
				getRadiusCode = function(source)
					-- TODO compute dr's component in each dx^i's, and scale our random number by that
					-- add some noise
					return '.2 - .003 * cos(1.+2.*M_PI*cos(1.+2.*M_PI*cos(1.+2.*M_PI*x.y/x.x)))'
				end,
				sources={
					{
						center={0, 0, 0}, 
						radius = .2,
						-- srhd solver requires the max velocity not to exceed 1 ...
						inside = inside,
					},
				},
			}(solver)
		end,
	},

	{
		name = 'self-gravitation test 2',
		initState = SelfGravProblem{ 
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
		initState = SelfGravProblem{ 
			sources = {
				{
					center = {-.25, 0, 0},
					radius = .1,
					inside = [[
						v.x = -.5 * x.y;
						v.y = .5 * x.x;
						rho = 1;
						P = 1;
					]],
				},
				{
					center = {.25, 0, 0},
					radius = .1,
					inside = [[
						v.x = -.5 * x.y;
						v.y = .5 * x.x;
						rho = 1;
						P = 1;
					]],
				},
			},
		},
	},
	
	{
		name = 'self-gravitation test 4',
		initState =  SelfGravProblem{
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
		initState = function(self, solver)
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
		initState = function(self, solver)
			return [[
	E = _real3(1,0,0);
	B = _real3(0, 1, lhs ? 1 : -1);
]]
		end,
	},

	{
		name = 'Maxwell scattering around cylinder',
		initState = function(self, solver)
			addMaxwellOscillatingBoundary(solver)
			return [[
	real3 xc = coordMap(x);
	if (real3_lenSq(xc) < .2*.2) {
		//conductivity = 0;
		permittivity = 10.;
	}
]]
		end,
	},

	{
		name = 'Maxwell wire',
		initState = function(self, solver)
			addMaxwellOscillatingBoundary(solver)
			
			local c = 299792458
			local s_in_m = 1 / c
			local G = 6.6740831e-11
			local kg_in_m = G / c^2
			local ke = 8.9875517873681764e+9
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
	conductivity = <?=clnumber(1/resistivities.air)?>;
	if (x.y >= -.1 && x.y < .1) {
		//if (x.x < mins.x * .9 + maxs.x * .1) 
			E.x = 1;
		conductivity = <?=clnumber(1/resistivities.copper)?>;
	}
]], 		{
				clnumber = clnumber,
				resistivities = resistivities,
			})
		end,
	},

	{
		name = 'two-fluid EMHD soliton ion',
		initState = function(self, solver)
			return [[
	const real L = 12.;
	rho = 1. + (real)exp((real)-25. * (real)fabs(x.x - L / (real)3.));
]]
		end,
	},
	{
		name = 'two-fluid EMHD soliton electron',
		initState = function(self, solver)
			return [[
	const real L = 12.;
	rho = 5. * (1. + (real)exp((real)-25. * (real)fabs(x.x - L / (real)3.)));
]]
		end,
	},
	{
		name = 'two-fluid EMHD soliton maxwell',
		initState = function(self, solver)
-- TODO			
			return [[

]]
		end,
	},
}:map(function(cl)
	return class(InitState, cl)
end)
return initStates

local function quadrantProblem(args)
	args.init = function(solver)
		solver.cfl[0] = .475
		for _,x in ipairs{'x', 'y', 'z'} do
			for _,minmax in ipairs{'min', 'max'} do
				solver.boundaryMethods[x..minmax][0] = solver.app.boundaryMethods:find'freeflow'-1
			end
		end
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

local function selfGravProblem(args)
	return function(solver)
		solver.useGravity = true
		solver.boundaryMethods.xmin[0] = solver.app.boundaryMethods:find'freeflow'-1
		solver.boundaryMethods.xmax[0] = solver.app.boundaryMethods:find'freeflow'-1
		solver.boundaryMethods.ymin[0] = solver.app.boundaryMethods:find'freeflow'-1
		solver.boundaryMethods.ymax[0] = solver.app.boundaryMethods:find'freeflow'-1
		solver.boundaryMethods.zmin[0] = solver.app.boundaryMethods:find'freeflow'-1
		solver.boundaryMethods.zmax[0] = solver.app.boundaryMethods:find'freeflow'-1

		return table{[[
	real distSq;
	rho = .1;
	P = 1;
	//v[i] = .1 * noise * crand();
]],	
			table.map(args.sources, function(source,i)
				local center = '(real4)('
					..range(3):map(function(i)
						return clnumber(source.center[i])
					end):concat', '..')'
				return ([[
	distSq = dot(x - center, x - center);
	if (dist < radius * radius) rho = 1; 
]]):gsub('center', center)
		:gsub('radius', tostring(source.radius))
			end),
		}:concat'\n'
	end
end

local initStates = {
	{
		name = 'Sod',
		init = function(solver)
			solver.eqn.guiVarsForName.heatCapacityRatio.value[0] = 7/5
			return [[
	rho = lhs ? 1 : .125;
	P = lhs ? 1 : .1;
]]
		end,
	},
	{
		name = 'Sedov',
		init = function(solver)
			return [[
	rho = 1;
	P = (i.x == gridSize.x/2 && i.y == gridSize.y/2 && i.z == gridSize.z/2) ? 1e+3 : 1;
]]
		end,
	},
	{
		name = 'constant',
		init = function(solver) 
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
		init = function(solver)
			return '	rho=2+x.x; P=1;'
		end,
	},
	{
		name = 'gaussian',
		init = function(solver)
			return [[
	real sigma = 1. / sqrt(10.);
	real xSq = dot(x,x);
	rho = exp(-xSq / (sigma*sigma)) + .1;
	P = 1 + .1 * (exp(-xSq / (sigma*sigma)) + 1) / ((heatCapacityRatio - 1.) * rho);
]]
		end,
	},
	{
		name = 'advect wave',
		init = function(solver)
			return [[
	real rSq = dot(x,x);
	rho = exp(-100*rSq) + 1.;
	v.x = 1;
	P = 1;
]]
		end,
	},
	-- http://www.cfd-online.com/Wiki/Explosion_test_in_2-D
	{
		name = 'sphere',
		init = function(solver)
			return [[
	real rSq = dot(x,x);
	bool inside = rSq < .2*.2;
	rho = inside ? 1 : .1;
	P = inside ? 1 : .1;
]]
		end,
	},
	{
		name = 'rarefaction wave',
		init = function(solver)
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
		init = function(solver)
			solver.eqn.guiVarsForName.heatCapacityRatio.value[0] = 4/3
			return [[
	rho = 1;
	v.x = 1. - 1e-5;
	P = (heatCapacityRatio - 1.) * rho * (1e-7 / sqrt(1. - v.x * v.x));
]]
		end,
	},
	{
		name = 'relativistic blast wave test problem 1',
		init = function(solver)
			solver.eqn.guiVarsForName.heatCapacityRatio.value[0] = 5/3
			return [[
	rho = lhs ? 10 : 1;
	P = (heatCapacityRatio - 1.) * rho * (lhs ? 2 : 1e-6);
]]
		end,
	},
	{
		name = 'relativistic blast wave test problem 2',
		init = function(solver)
			solver.eqn.guiVarsForName.heatCapacityRatio.value[0] = 5/3
			return [[
	rho = 1;
	P = lhs ? 1000 : .01;
]]
		end,
	},
	{
		name = 'relativistic blast wave interaction',
		init = function(solver)
			return [[
	real xL = .9 * mins_x + .1 * maxs_x;
	real xR = .1 * mins_x + .9 * maxs_x;
	rho = 1;
	P = x.x < xL ? 1000 : (x.x > xR ? 100 : .01);
]]
		end,
	},

	{
		name = 'Colella-Woodward',
		init = function(solver)
			for _,x in ipairs{'x', 'y', 'z'} do
				for _,minmax in ipairs{'min', 'max'} do
					solver.boundaryMethods[x..minmax][0] = solver.app.boundaryMethods:find'freeflow'-1
				end
			end
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
		init = function(solver)
			for _,x in ipairs{'x', 'y', 'z'} do
				for _,minmax in ipairs{'min', 'max'} do
					solver.boundaryMethods[x..minmax][0] = solver.app.boundaryMethods:find'periodic'-1
				end
			end
			return [[
	bool inside = (x.y > -.25 && x.y < .25);
	real theta = (x.x - mins.x) / (maxs.x - mins.x) * 2. * M_PI;
#if dim == 3
	theta *= (x.z - mins.z) / (maxs.z - mins.z);
#endif
	real noise = (maxs.x - mins.x) * 1e-4;
	rho = inside ? 2 : 1;
	v.x = cos(theta) * noise + (inside ? -.5 : .5);
#if dim == 2
	v.y = sin(theta) * noise;
#endif
#if dim == 3
	v.z = sin(theta) * noise;
#endif
	P = 2.5;
]]
		end,
	},
	
	--http://www.astro.virginia.edu/VITA/ATHENA/dmr.html
	{
		name = 'double mach reflection',
		init = function(solver)
			-- I am not correctly modeling the top boundary
			solver.mins = vec3(0,0,0)
			solver.maxs = vec3(4,1,1)
			solver.boundaryMethods.xmin[0] = solver.app.boundaryMethods:find'freeflow'-1
			solver.boundaryMethods.xmax[0] = solver.app.boundaryMethods:find'freeflow'-1
			solver.boundaryMethods.ymin[0] = solver.app.boundaryMethods:find'mirror'-1
			solver.boundaryMethods.ymax[0] = solver.app.boundaryMethods:find'freeflow'-1
			solver.boundaryMethods.zmin[0] = solver.app.boundaryMethods:find'mirror'-1
			solver.boundaryMethods.zmax[0] = solver.app.boundaryMethods:find'mirror'-1
			solver.guiVarsForName.heatCapacityRatio.value[0] = 7/5
			return table{
	'#define sqrt1_3 '..clnumber(math.sqrt(1/3)),
	[[
	bool inside = x.x < x.y * sqrt1_3;
	if (inside) {
		rho = 8;
		P = 116.5;
		v.x = 8.25 * cos(30. * M_PI / 180.),
		v.y = -8.25 * sin(30. * M_PI / 180.),
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
		init = function(solver)
			solver.boundaryMethods.xmin[0] = solver.app.boundaryMethods:find'mirror'-1
			solver.boundaryMethods.xmax[0] = solver.app.boundaryMethods:find'mirror'-1
			solver.boundaryMethods.ymin[0] = solver.app.boundaryMethods:find'mirror'-1
			solver.boundaryMethods.ymax[0] = solver.app.boundaryMethods:find'freeflow'-1	-- TODO none
			solver.boundaryMethods.zmin[0] = solver.app.boundaryMethods:find'mirror'-1
			solver.boundaryMethods.zmax[0] = solver.app.boundaryMethods:find'mirror'-1
			return [[
	rho = 1;
	v.x = x.y > .45 ? 1 : 0;
	P = 1;
]]
		end,
	},

	{
		name='shock bubble interaction',
		init = function(solver)
			solver.boundaryMethods.xmin[0] = solver.app.boundaryMethods:find'freeflow'-1
			solver.boundaryMethods.xmax[0] = solver.app.boundaryMethods:find'freeflow'-1
			solver.boundaryMethods.ymin[0] = solver.app.boundaryMethods:find'freeflow'-1
			solver.boundaryMethods.ymax[0] = solver.app.boundaryMethods:find'freeflow'-1
			solver.boundaryMethods.zmin[0] = solver.app.boundaryMethods:find'freeflow'-1
			solver.boundaryMethods.zmax[0] = solver.app.boundaryMethods:find'freeflow'-1
			return [[
	const real waveX = .45;
	const real2 bubbleCenter = (real2)(0,0);
	real bubbleRadius = .2;
	real bubbleRSq = dot(x - bubbleCenter, x - bubbleCenter);
	rho = x.x < waveX ? 1. : (bubbleRSq < bubbleRadius*bubbleRadius ? .1 : 1);
	P = x.x < waveX ? 1 : .1;
	v.x = x.x < waveX ? 0 : -.5;
]]
		end,
	},

	-- gravity potential test - equilibrium - Rayleigh-Taylor (still has an shock wave ... need to fix initial conditions?)

	{
		name = 'self-gravitation test 1',
		init = selfGravProblem{
			solver = solver,
			sources={
				{center={0, 0, 0}, radius = .2},
			},
		},
	},

	{
		name = 'self-gravitation test 1 spinning',
		init = selfGravProblem{
			sources={
				{
					center={0, 0, 0}, 
					radius = .2,
					inside = function(dx,dy,dz)
						return buildStateEuler{
							x=x, y=y, z=z,
							velocityX = -10 * dy,
							velocityY = 10 * dx,
							pressure = 1,
							density = 1,
						}
					end},
			},
		},
	},

	{
		name = 'self-gravitation test 2',
		init = selfGravProblem{ 
			sources={
				{
					center = {-.25, 0, 0},
					radius = .1,
					inside = function(dx,dy,dz)
						return buildStateEuler{
							x=x, y=y, z=z,
							pressure = 1,
							density = 1,
						}
					end,
				},
				{
					center = {.25, 0, 0},
					radius = .1,
					inside = function(dx,dy,dz)
						return buildStateEuler{
							x=x, y=y, z=z,
							pressure = 1,
							density = 1,
						}
					end,
				},
			},
		},
		--[[ TODO
		mx = -5 * rho * y
		my = 5 * rho * x
		return rho,mx,my,mz,eTotal,bx,by,bz
		--]]
	},

	{
		name = 'self-gravitation test 4',
		init =  selfGravProblem{
			sources={
				{center={.25, .25, 0}, radius = .1},
				{center={-.25, .25, 0}, radius = .1},
				{center={.25, -.25, 0}, radius = .1},
				{center={-.25, -.25, 0}, radius = .1},
			},
		},
	},
}

return initStates

local table = require 'ext.table'
local range = require 'ext.range'
local clnumber = require 'cl.obj.number'
local template = require 'template'

local function quadrantProblem(args)
	args.init = function(solver)
		solver.cfl = .475
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

		return template([[
	rho = .1;
	P = 1;
	//v[i] = .1 * noise * crand();
	<? for i,source in ipairs(sources) do ?>{
		real3 delta = real3_sub(x, _real3(
			<?=clnumber(source.center[1])?>,
			<?=clnumber(source.center[2])?>,
			<?=clnumber(source.center[3])?>));
		real distSq = real3_lenSq(delta);
		if (distSq < <?=clnumber(source.radius * source.radius)?>) rho = 1; 
	}<? end ?>
]], 	{
			sources = args.sources,
			clnumber = clnumber,
		})
	end
end

local initStates = {
	{
		name = 'constant',
		init = function(solver)
			if solver.eqn.guiVarsForName.heatCapacityRatio then	
				solver.eqn.guiVarsForName.heatCapacityRatio.value[0] = 7/5
			end
			return [[
	rho = 1;
	P = 1;
]]
		end,
	},
	{
		name = 'Sod',
		init = function(solver)
			if solver.eqn.guiVarsForName.heatCapacityRatio then	
				solver.eqn.guiVarsForName.heatCapacityRatio.value[0] = 7/5
			end
			return [[
	rho = lhs ? 1. : .125;
	P = lhs ? 1. : .1;
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
	-- http://www.astro.uni-bonn.de/~jmackey/jmac/node7.html
	-- http://www.astro.princeton.edu/~jstone/Athena/tests/brio-wu/Brio-Wu.html
	{
		name = 'Brio-Wu',
		init = function(solver)
			if solver.eqn.guiVarsForName.heatCapacityRatio then	
				solver.eqn.guiVarsForName.heatCapacityRatio.value[0] = 2
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
		init = function(solver)
			if solver.eqn.guiVarsForName.heatCapacityRatio then	
				solver.eqn.guiVarsForName.heatCapacityRatio.value[0] = 5/3
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
			solver.cfl = .5	-- needs a slower cfl
			if solver.eqn.guiVarsForName.heatCapacityRatio then	
				solver.eqn.guiVarsForName.heatCapacityRatio.value[0] = 4/3
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
		init = function(solver)
			solver.cfl = .5	-- needs a slower cfl
			if solver.eqn.guiVarsForName.heatCapacityRatio then	
				solver.eqn.guiVarsForName.heatCapacityRatio.value[0] = 5/3
			end
			return [[
	rho = lhs ? 10 : 1;
	P = (heatCapacityRatio - 1.) * rho * (lhs ? 2 : 1e-6);
]]
		end,
	},
	{
		name = 'relativistic blast wave test problem 2',
		init = function(solver)
			solver.cfl = .5	-- needs a slower cfl
			if solver.eqn.guiVarsForName.heatCapacityRatio then	
				solver.eqn.guiVarsForName.heatCapacityRatio.value[0] = 5/3
			end
			return [[
	rho = 1;
	P = lhs ? 1000 : .01;
]]
		end,
	},
	{
		name = 'relativistic blast wave interaction',
		init = function(solver)
			solver.cfl = .5	-- needs a slower cfl
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
			for i,x in ipairs{'x', 'y', 'z'} do
				for _,minmax in ipairs{'min', 'max'} do
					solver.boundaryMethods[x..minmax][0] = 
						--i == solver.dim and solver.app.boundaryMethods:find'freeflow'-1 or 
						solver.app.boundaryMethods:find'periodic'-1
				end
			end
			
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
	{
		name = 'Rayleigh-Taylor',
		init = function()
			local xs = {'x', 'y', 'z'}
			
			local xmid = (solver.mins + solver.maxs) * .5
			
			-- triple the length along the interface dimension
			local k = xs[solver.dim]
			solver.mins[k] = xmid[k] + (solver.mins[k] - xmid[k]) * 3
			solver.maxs[k] = xmid[k] + (solver.maxs[k] - xmid[k]) * 3

			-- triple resolution along interface dimension
			--size[k] = size[k] * 3
			-- (can it handle npo2 sizes?)
		
			for i,x in ipairs(xs) do
				for _,minmax in ipairs{'min', 'max'} do
					solver.boundaryMethods[x..minmax][0] = 
						i == solver.dim 
						and solver.app.boundaryMethods:find'mirror'-1
						or solver.app.boundaryMethods:find'periodic'-1
				end
			end
		
			-- TODO incorporate this into the Euler model ..
			local externalForce = {0, 1, 0}
			
			initState = function(x,y,z)
				local xs = {x,y,z}
				local top = xs[k] > xmid[k]
				local potentialEnergy = 0	-- minPotentialEnergy
				for k=1,#size do
					potentialEnergy = potentialEnergy + (xs[k] - xmin[k]) * externalForce[k]
				end
				local density = top and 2 or 1
				return buildStateEuler{
					noise = .001,
					density = density,
					potentialEnergy = potentialEnergy,
					pressure = 2.5 - density * potentialEnergy,
					-- or maybe it is ... pressure = (gamma - 1) * density * (2.5 - potentialEnergy)
				}
			end
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
			if solver.eqn.guiVarsForName.heatCapacityRatio then	
				solver.guiVarsForName.heatCapacityRatio.value[0] = 7/5
			end
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

	{
		name = 'Maxwell default',
		init = function(solver)
			return [[
	const real L = 12.;
	E = _real3(1,0,0);
	B = _real3(0, 1, lhs ? 1 : -1);
]]
		end,
	},

	{
		name = 'two-fluid EMHD soliton ion',
		init = function(solver)
			return [[
	const real L = 12.;
	rho = 1. + (real)exp((real)-25. * (real)fabs(x.x - L / (real)3.));
]]
		end,
	},
	{
		name = 'two-fluid EMHD soliton electron',
		init = function(solver)
			return [[
	const real L = 12.;
	rho = 5. * (1. + (real)exp((real)-25. * (real)fabs(x.x - L / (real)3.)));
]]
		end,
	},
	{
		name = 'two-fluid EMHD soliton maxwell',
		init = function(solver)
-- TODO			
			return [[

]]
		end,
	},
}

return initStates

local class = require 'ext.class'
local table = require 'ext.table'
local InitCond = require 'hydro.init.init'

local NLSInitCond = class(InitCond)

local initConds = table{
	{
		name = 'Gaussian',
		mins = {.1, .1, .1},
		maxs = {4, 4, 4},
		guiVars = {
			{name='A', value=10},
		},
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			-- TODO custom boundary.  rhs is set to zero.  lhs is U[-2] = U[2], U[-1] = U[1], and U[0] is not modified
			solver:setBoundaryMethods'freeflow'
			return [[
	q = cplx_from_real(initCond->A * exp(-r * r));
]]
		end,
	},
	{
		name = 'Ring',
		mins = {.1,.1,.1},
		maxs = {4,4,4},
		guiVars = {
			{name='A', value=8},
		},
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			solver:setBoundaryMethods'freeflow'
			return [[
	q = cplx_from_real(initCond->A * r * r * exp(-r * r));
]]
		end,
	},
	{
		name = 'Oscillatory',
		mins = {.1, .1, .1},
		maxs = {4, 4, 4},
		guiVars = {	
			{name='A', value=4},
			{name='alpha', value=10},
		},
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			solver:setBoundaryMethods'freeflow'
			return [[
	real magn = initCond->A * exp(-r * r);
	real theta = -initCond->alpha * r * r;
	q = cplx_from_real(
		cos(theta) * magn,
		sin(theta) * magn
	);
]]
		end,
	},

	----------------------------------------------------------
	---- piggybacking for my wave-finite-difference here: ----
	----------------------------------------------------------
	
	{
		name = 'Wave-FD Gaussian',
		mins = {.3, .3, .3}, 
		maxs = {20.3, 20.3, 20.3},-- paper uses .3 + 20	
		--maxs = {5.3, 5.3, 5.3},	
		-- from 2014 Oliveira et al PRD
		-- paper says grid uses rmax=200
		-- paper also says r = rmin + j h  "for h the resolution of the grid"
		-- paper also says h = 1/30 ... 1/2000
		-- so is rmax fixed and h determined by (rmin-rmax)/n, 
		--  or is h fixed and is rmax = rmin + n h ?
		guiVars = {
			{name='r0', value=2},
			{name='sigma', value=.25},
		},
		getInitCondCode = function(self)
			return [[
	real const rmin = solver->mins.x;
	real const drmin = r - rmin;
	real const dr0_over_sigma = (r - initCond->r0) / initCond->sigma;
	q = cplx_from_real(drmin * drmin * exp(-.5 * dr0_over_sigma * dr0_over_sigma));
]]
		end,
	},

	{
		name = 'Wave-FD Bessel',
		mins = {.3, .3, .3}, 
		maxs = {20.3, 20.3, 20.3},	
		depends = {'Bessel'},
		getInitCondCode = function(self)
			return [[
	//q = cplx_from_real(BESSJ0(x.x));
	// bessel deriv
	q = cplx_from_real(BESSJ1(x.x));
]]
		end,
	},

	-- made for schrodinger-fd
	-- https://colab.research.google.com/drive/10OcwtxjXb9xiscx-6eQiG1hoUFYLLqvq?usp=sharing#scrollTo=pp8YgtAZ3onC
	{
		name = 'electron in potential',
		mins = {0, 0, 0},
		maxs = {.01, .01, .01},
		getInitCondCode = function(self)
print'TODO fix the numerical precision errors'			
			local solver = assert(self.solver)
			return solver.eqn:template[[
//// MODULE_DEPENDS: sqr <?=coordMap?>
	real const xmax = (real)solver->gridSize.x;
	real const ymax = (real)solver->gridSize.y;
	real3 const xc = coordMap(x);
	real const A = exp(-(( sqr(xc.y-ymax/2) ))/1000.) * exp(-(( sqr(xc.x-xmax/2-xmax/4-xmax/8) ))/1000);
	q.re = sin((60.*M_PI*xc.x)/xmax);
	q.im = -cos((60.*M_PI*xc.x)/xmax+M_PI);

    V = cplx_from_real(exp(-( sqr(xc.x-350)/200. + sqr(xc.y-ymax/2)/200.)) * 4e-29);

]]
		end,
	},
}:mapi(function(cl)
	return class(NLSInitCond, cl)
end)

function NLSInitCond:getList()
	return initConds
end

return NLSInitCond 

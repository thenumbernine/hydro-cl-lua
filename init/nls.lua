local class = require 'ext.class'
local table = require 'ext.table'
local InitCond = require 'init.init'

return table{
	{
		name = 'Gaussian',
		mins = {.1, .1, .1},
		maxs = {4, 4, 4},
		guiVars = {
			{name='nls_A', value=10},
		},
		initState = function(self, solver)
			-- TODO custom boundary.  rhs is set to zero.  lhs is U[-2] = U[2], U[-1] = U[1], and U[0] is not modified
			solver:setBoundaryMethods'freeflow'
			return [[
	q = cplx_from_real(solver->nls_A * exp(-r * r));
]]
		end,
	},
	{
		name = 'Ring',
		mins = {.1,.1,.1},
		maxs = {4,4,4},
		guiVars = {
			{name='nls_A', value=8},
		},
		initState = function(self, solver)
			solver:setBoundaryMethods'freeflow'
			return [[
	q = cplx_from_real(solver->nls_A * r * r * exp(-r * r));
]]
		end,
	},
	{
		name = 'Oscillatory',
		mins = {.1, .1, .1},
		maxs = {4, 4, 4},
		guiVars = {	
			{name='nls_A', value=4},
			{name='nls_alpha', value=10},
		},
		initState = function(self, solver)
			solver:setBoundaryMethods'freeflow'
			return [[
	real magn = solver->nls_A * exp(-r * r);
	real theta = -nls_alpha * r * r;
	q = cplx_from_real(
		cos(theta) * magn,
		sin(theta) * magn
	);
]]
		end,
	},

	-- piggybacking for my wave-finite-difference here:
	{
		name = 'Wave-FD Gaussian',
		mins = {.3, .3, .3}, 
		
		maxs = {20.3, 20.3, 20.3},	
		-- paper says grid uses rmax=200
		-- paper also says r = rmin + j h  "for h the resolution of the grid"
		-- paper also says h = 1/30 ... 1/2000
		-- so is rmax fixed and h determined by (rmin-rmax)/n, 
		--  or is h fixed and is rmax = rmin + n h ?
		guiVars = {
			-- TODO put these in eqn/wave-fd.lua
			{name='m', value=2},
			{name='C', value=.5},
		
			{name='init_r0', value=2},
			{name='init_sigma', value=.25},
		},
		initState = function(self, solver)
			return [[
	real rmin = solver->mins.x;
	real drmin = r - rmin;
	real dr0_over_sigma = (r - solver->init_r0) / solver->init_sigma;
	q = cplx_from_real(drmin * drmin * exp(-.5 * dr0_over_sigma * dr0_over_sigma));
]]
		end,
	},
}:map(function(cl)
	return class(InitCond, cl)
end)

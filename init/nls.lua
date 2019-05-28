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
	re = solver->nls_A * exp(-r * r);
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
	re = solver->nls_A * r * r * exp(-r * r);
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
	re = cos(theta) * magn;
	im = sin(theta) * magn;
]]
		end,
	},
}:map(function(cl)
	return class(InitCond, cl)
end)

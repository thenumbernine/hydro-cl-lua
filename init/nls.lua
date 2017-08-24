local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local clnumber = require 'cl.obj.number'
local vec3 = require 'vec.vec3'
local InitCond = require 'init.init'

return table{
	{
		name = 'Gaussian',
		init = function(self, solver)
			solver.eqn:addGuiVar{name='nls_A', value=10}
			solver.mins = vec3(.1,.1,.1)
			solver.maxs = vec3(4,4,4)
			-- TODO custom boundary.  rhs is set to zero.  lhs is U[-2] = U[2], U[-1] = U[1], and U[0] is not modified
			solver:setBoundaryMethods'freeflow'
		end,
		initState = function(self, solver)
			return [[
	re = nls_A * exp(-r * r);
]]
		end,
	},
	{
		name = 'Ring',
		init = function(self, solver)
			solver.eqn:addGuiVar{name='nls_A', value=8}
			solver.mins = vec3(.1,.1,.1)
			solver.maxs = vec3(4,4,4)
			solver:setBoundaryMethods'freeflow'
		end,
		initState = function(self, solver)
			return [[
	re = nls_A * r * r * exp(-r * r);
]]
		end,
	},
	{
		name = 'Oscillatory',
		init = function(self, solver)
			solver.eqn:addGuiVar{name='nls_A', value=4}
			solver.eqn:addGuiVar{name='nls_alpha', value=10}
			solver.mins = vec3(.1,.1,.1)
			solver.maxs = vec3(4,4,4)
			solver:setBoundaryMethods'freeflow'
		end,
		initState = function(self, solver)
			return [[
	real magn = nls_A * exp(-r * r);
	real theta = -nls_alpha * r * r;
	re = cos(theta) * magn;
	im = sin(theta) * magn;
]]
		end,
	},
}:map(function(cl)
	return class(InitCond, cl)
end)

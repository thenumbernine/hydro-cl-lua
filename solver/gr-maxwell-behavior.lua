local GRBehavior = require 'solver.gr-behavior'
local GRMaxwellEqn = require 'eqn.gr-maxwell'
return function(parent)
	local cl = GRBehavior(parent)
	cl.eqnName = 'gr-maxwell'
	return cl
end

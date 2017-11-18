local GRBehavior = require 'solver.gr-behavior'
local GRMaxwellEqn = require 'eqn.gr-maxwell'
return function(parent)
	return GRBehavior(parent, GRMaxwellEqn)
end

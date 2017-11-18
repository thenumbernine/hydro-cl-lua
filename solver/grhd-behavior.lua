local GRBehavior = require 'solver.gr-behavior'
local RHDBehavior = require 'solver.rhd-behavior'
local GRHDEqn = require 'eqn.grhd'
return function(parent)
	return GRBehavior(RHDBehavior(parent), GRHDEqn)
end

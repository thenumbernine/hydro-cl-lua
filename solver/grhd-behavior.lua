local class = require 'ext.class'
local GRBehavior = require 'solver.gr-behavior'
local RHDBehavior = require 'solver.rhd-behavior'
return function(parent)
	local template = class(GRBehavior(RHDBehavior(parent)))
	template.eqnName = 'grhd'
	return template
end

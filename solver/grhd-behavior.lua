local class = require 'ext.class'
local GRBehavior = require 'solver.gr-behavior'
return function(parent)
	local template = class(GRBehavior(parent))
	template.eqnName = 'grhd'
	return template
end

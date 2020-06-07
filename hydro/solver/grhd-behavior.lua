local class = require 'ext.class'
local GRBehavior = require 'hydro.solver.gr-behavior'
return function(parent)
	local template = class(GRBehavior(parent))
	template.eqnName = 'grhd'
	return template
end

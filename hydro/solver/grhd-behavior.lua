local GRBehavior = require 'hydro.solver.gr-behavior'
return function(parent)
	local template = GRBehavior(parent):subclass()
	template.eqnName = 'grhd'
	return template
end

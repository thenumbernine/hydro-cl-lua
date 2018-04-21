local class = require 'ext.class'
local RHDBehavior = require 'solver/rhd-behavior'
return function(parent)
	local template = class(RHDBehavior(parent))
	template.eqnName = 'srhd'
	return template
end

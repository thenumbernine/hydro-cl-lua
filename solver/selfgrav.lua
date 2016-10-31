local class = require 'ext.class'
local SelfGravitationBehavior = function(parent)
	local template = class(parent)

	function template:update(...)
		template.super.update(self, ...)

		-- TODO calculate potential from density here
		-- TODO integrate potential gradient to velocity here 
	end

	return template
end

return SelfGravitationBehavior 

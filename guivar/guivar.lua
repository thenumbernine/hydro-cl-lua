local class = require 'ext.class'

local GuiVar = class()

function GuiVar:init(args)
	self.name = assert(args.name)
end

return GuiVar

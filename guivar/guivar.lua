local class = require 'ext.class'

local GuiVar = class()

function GuiVar:init(args)
	self.name = assert(args.name)
end

function GuiVar:refresh(value, solver)
	print('refreshing '..self.name..' = '..value)
	solver:refreshCodePrefix()
end

return GuiVar

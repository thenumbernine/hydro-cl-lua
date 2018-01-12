local class = require 'ext.class'

local GuiVar = class()

function GuiVar:init(args)
	self.name = assert(args.name)
	self.onChange = args.onChange
end

function GuiVar:refresh(value, solver)
	print('refreshing '..self.name..' = '..tostring(value))
	if self.onChange then self:onChange(value, solver) end
	solver:refreshCodePrefix()
end

return GuiVar

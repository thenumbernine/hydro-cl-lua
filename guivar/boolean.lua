local ig = require 'ffi.imgui'
local tooltip = require 'tooltip'
local class = require 'ext.class'
local GuiVar = require 'guivar.guivar'

local GuiBoolean = class(GuiVar)

GuiBoolean.ctype = 'int'

function GuiBoolean:init(args)
	GuiBoolean.super.init(self, args)
	self.value = not not args.value
end

function GuiBoolean:updateGUI(solver)
	if tooltip.checkboxTable(self.name, self, 'value') then
		self:refresh(self.value, solver)
	end
end

-- compile-time
function GuiBoolean:getCode()
	return '#define '..self.name..' '..(self.value and 1 or 0)
end

return GuiBoolean

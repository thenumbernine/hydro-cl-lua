local ig = require 'imgui'
local GuiVar = require 'hydro.guivar.guivar'

local GuiBoolean = GuiVar:subclass()

GuiBoolean.ctype = 'int'

function GuiBoolean:init(args)
	GuiBoolean.super.init(self, args)
	self.value = not not args.value
end

function GuiBoolean:updateGUI(solver)
	if ig.luatableTooltipCheckbox(self.name, self, 'value') then
		self:refresh(solver)
	end
end

-- compile-time
function GuiBoolean:getCode()
	return '#define '..self.name..' '..(self.value and 1 or 0)
end

return GuiBoolean

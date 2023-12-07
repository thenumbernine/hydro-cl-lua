local ig = require 'imgui'
local clnumber = require 'cl.obj.number'
local GuiVar = require 'hydro.guivar.guivar'

local GuiNumber = GuiVar:subclass()

GuiNumber.ctype = 'real' 

function GuiNumber:init(args)
	GuiNumber.super.init(self, args)
	self.value = args.value or 0
end

function GuiNumber:updateGUI(solver)
	if ig.luatableTooltipInputFloatAsText(self.name, self, 'value', ig.ImGuiInputTextFlags_EnterReturnsTrue) then
		self:refresh(solver)
	end
end

-- compile-time
function GuiNumber:getCode()
	return '#define '..self.name..' '..clnumber(self.value)
end

return GuiNumber

local ffi = require 'ffi'
local ig = require 'imgui'
local GuiVar = require 'hydro.guivar.guivar'

local GuiInt = GuiVar:subclass()

GuiInt.ctype = 'int' 

function GuiInt:init(args)
	GuiInt.super.init(self, args)
	self.value = args.value or 0
end

function GuiInt:updateGUI(solver)
	if ig.luatableTooltipInputInt(self.name, self, 'value', 1, 100, ig.ImGuiInputTextFlags_EnterReturnsTrue) then
		self:refresh(solver)
	end
end

-- compile-time
function GuiInt:getCode()
	return '#define '..self.name..' '..self.value
end

return GuiInt

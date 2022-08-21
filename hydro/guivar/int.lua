local ffi = require 'ffi'
local ig = require 'imgui'
local class = require 'ext.class'
local GuiVar = require 'hydro.guivar.guivar'

local GuiInt = class(GuiVar)

GuiInt.ctype = 'int' 

function GuiInt:init(args)
	GuiInt.super.init(self, args)
	self.value = args.value or 0
end

function GuiInt:updateGUI(solver)
	if ig.luatableTooltipInputInt(self.name, self, 'value', 1, 100, ig.ImGuiInputTextFlags_EnterReturnsTrue) then
		self:refresh(self.value, solver)
	end
end

-- compile-time
function GuiInt:getCode()
	return '#define '..self.name..' '..self.value
end

return GuiInt

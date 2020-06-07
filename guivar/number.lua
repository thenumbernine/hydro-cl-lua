local ig = require 'ffi.imgui'
local tooltip = require 'hydro.tooltip'
local class = require 'ext.class'
local GuiVar = require 'guivar.guivar'
local clnumber = require 'cl.obj.number'

local GuiNumber = class(GuiVar)

GuiNumber.ctype = 'real' 

function GuiNumber:init(args)
	GuiNumber.super.init(self, args)
	self.value = args.value or 0
end

function GuiNumber:updateGUI(solver)
	if tooltip.numberTable(self.name, self, 'value', ig.ImGuiInputTextFlags_EnterReturnsTrue) then
		self:refresh(self.value, solver)
	end
end

-- compile-time
function GuiNumber:getCode()
	return '#define '..self.name..' '..clnumber(self.value)
end

return GuiNumber

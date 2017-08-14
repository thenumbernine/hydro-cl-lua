local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local tooltip = require 'tooltip'
local class = require 'ext.class'
local GuiVar = require 'guivar.guivar'
local clnumber = require 'cl.obj.number'

local GuiFloat = class(GuiVar)

function GuiFloat:init(args)
	GuiFloat.super.init(self, args)
	self.value = args.value or 0
end

function GuiFloat:getCode()
	return '#define '..self.name..' '..clnumber(self.value)
end

function GuiFloat:updateGUI(solver)
	if tooltip.numberTable(self.name, self, 'value', ig.ImGuiInputTextFlags_EnterReturnsTrue) then
		self:refresh(self.value, solver)
	end
end

return GuiFloat

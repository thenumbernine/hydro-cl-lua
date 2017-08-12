local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local tooltip = require 'tooltip'
local class = require 'ext.class'
local GuiVar = require 'guivar.guivar'

local GuiInt = class(GuiVar)

function GuiInt:init(args)
	GuiInt.super.init(self, args)
	self.value = args.value or 0
end

function GuiInt:getCode()
	return '#define '..self.name..' '..self.value
end

function GuiInt:updateGUI(solver)
	if tooltip.intTable(self.name, self, 'value', 1, 100, ig.ImGuiInputTextFlags_EnterReturnsTrue) then
		self:refresh(self.value, solver)
	end
end

return GuiInt

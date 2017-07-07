local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local tooltip = require 'tooltip'
local class = require 'ext.class'
local GuiVar = require 'guivar.guivar'

local GuiCombo = class(GuiVar)

function GuiCombo:init(args)
	GuiCombo.super.init(self, args)
	self.value = args.value or 1
	self.options = assert(args.options)
end

function GuiCombo:getCode()
	return '#define '..self.name..' '..self.options[self.value]
end

function GuiCombo:updateGUI(solver)
	if tooltip.comboTable(self.name, self, 'value', self.options) then
		print('refreshing '..self.name..' = '..self.options[self.value])
		solver:refreshSolverProgram()
		solver:refreshDisplayProgram()
	end
end

return GuiCombo

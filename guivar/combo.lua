local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local GuiVar = require 'guivar.guivar'

local GuiCombo = class(GuiVar)

function GuiCombo:init(args)
	GuiCombo.super.init(self, args)
	self.value = ffi.new('int[1]', args.value or 0)
	self.options = assert(args.options)
end

function GuiCombo:getCode()
	return '#define '..self.name..' '..self.options[self.value[0]+1]
end

function GuiCombo:updateGUI(solver)
	if ig.igCombo(self.name, self.value, self.options) then
		print('refreshing '..self.name..' = '..self.options[self.value[0]+1])
		solver:refreshSolverProgram()
		solver:refreshDisplayProgram()
	end
end

return GuiCombo

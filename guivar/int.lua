local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local GuiVar = require 'guivar.guivar'

local GuiInt = class(GuiVar)

function GuiInt:init(args)
	GuiInt.super.init(self, args)
	self.value = ffi.new('int[1]', args.value or 0)
end

function GuiInt:getCode()
	return '#define '..self.name..' '..self.value[0]
end

function GuiInt:updateGUI(solver)
	if ig.igInputInt(self.name, self.value, 1, 100, ig.ImGuiInputTextFlags_EnterReturnsTrue) then
		print('refreshing '..self.name..' = '..self.value[0])
		solver:refreshSolverProgram()
		solver:refreshDisplayProgram()
	end
end

return GuiInt

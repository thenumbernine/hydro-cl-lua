local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local GuiVar = require 'guivar.guivar'
local clnumber = require 'clnumber'

local GuiFloat = class(GuiVar)

function GuiFloat:init(args)
	GuiFloat.super.init(self, args)
	self.value = ffi.new('float[1]', args.value or 0)
end

function GuiFloat:getCode()
	return '#define '..self.name..' '..clnumber(self.value[0])
end

function GuiFloat:updateGUI(solver)
	if ig.igInputFloat(self.name, self.value, 0, 0, -1, ig.ImGuiInputTextFlags_EnterReturnsTrue) then
		print('refreshing '..self.name..' = '..self.value[0])
		solver:refreshSolverProgram()
		solver:refreshDisplayProgram()
	end
end

return GuiFloat

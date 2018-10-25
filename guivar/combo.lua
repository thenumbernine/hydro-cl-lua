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

function GuiCombo:updateGUI(solver)
	if tooltip.comboTable(self.name, self, 'value', self.options) then
		self:refresh(self.options[self.value], solver)
	end
end

-- TODO how about another function for generating the enum code, whether or not this is a compile-time or runtime variable

-- compile-time
function GuiCombo:getCode()
	return '#define '..self.name..' '..self.options[self.value]
end

-- run-time
function GuiCombo:addToSolver(solver)
	solver.solverVars:insert{[self.name] = 'int'}
end

function GuiCombo:setToSolver(solver)
	solver.solverPtr[self.name] = self.value
end

return GuiCombo

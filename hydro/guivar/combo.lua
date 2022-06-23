local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local class = require 'ext.class'
local tooltip = require 'hydro.tooltip'
local GuiVar = require 'hydro.guivar.guivar'

local GuiCombo = class(GuiVar)

GuiCombo.ctype = 'int'

function GuiCombo:init(args)
	GuiCombo.super.init(self, args)
	self.value = args.value or 1
	self.options = assert(args.options)
end

function GuiCombo:updateGUI(solver)
	if tooltip.comboTable(self.name, self, 'value', self.options) then
		self:refresh(self:getValue(), solver)
	end
end

-- TODO how about another function for generating the enum code, whether or not this is a compile-time or runtime variable

-- compile-time
function GuiCombo:getCode()
	return '#define '..self.name..' '..self:getValue()
end

-- ok .value is really the key ...
-- ... and :getValue() is the value
function GuiCombo:getValue()
	return self.options[self.value]
end

return GuiCombo

local ffi = require 'ffi'
local ig = require 'imgui'
local GuiVar = require 'hydro.guivar.guivar'

local GuiCombo = GuiVar:subclass()

GuiCombo.ctype = 'int'

function GuiCombo:init(args)
	GuiCombo.super.init(self, args)
	self.value = args.value or 1
	self.options = assert(args.options)
end

function GuiCombo:updateGUI(solver)
	if ig.luatableTooltipCombo(self.name, self, 'value', self.options) then
		self:refresh(solver)
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

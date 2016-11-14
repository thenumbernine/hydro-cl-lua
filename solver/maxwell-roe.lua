local class = require 'ext.class'
local Roe = require 'solver.roe'
local NoDiv = require 'solver.nodiv'

local MaxwellRoe = class(NoDiv(Roe))

function MaxwellRoe:init(...)
	MaxwellRoe.super.init(self, ...)

	select(2, self.displayVars:find(nil, function(var) return var.name == 'U_div_B' end)).enabled[0] = true
	select(2, self.displayVars:find(nil, function(var) return var.name == 'U_div_B' end)).heatMapFixedRangePtr[0] = false 
	select(2, self.displayVars:find(nil, function(var) return var.name == 'U_div_E' end)).enabled[0] = true 
	select(2, self.displayVars:find(nil, function(var) return var.name == 'U_div_E' end)).heatMapFixedRangePtr[0] = false 
	select(2, self.displayVars:find(nil, function(var) return var.name == 'ePot_0' end)).enabled[0] = true 
	select(2, self.displayVars:find(nil, function(var) return var.name == 'ePot_0' end)).heatMapFixedRangePtr[0] = false 
	self:refreshDisplayProgram()
end

function MaxwellRoe:createEqn()
	self.eqn = require 'eqn.maxwell'(self)
end

return MaxwellRoe

local class = require 'ext.class'
local Roe = require 'solver.roe'
local NoDiv = require 'solver.nodiv'

local MaxwellRoe = class(NoDiv(Roe))

function MaxwellRoe:createEqn()
	self.eqn = require 'eqn.maxwell'(self)
end

return MaxwellRoe

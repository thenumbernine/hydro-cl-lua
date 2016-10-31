local class = require 'ext.class'
local Roe = require 'solver.roe'
local SelfGravitationBehavior = require 'solver.selfgrav'

--local EulerRoe = class(SelfGravitationBehavior(Roe))
local EulerRoe = class(Roe)

function EulerRoe:createEqn()
	self.eqn = require 'eqn.euler'(self)
end

return EulerRoe

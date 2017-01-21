local ffi = require 'ffi'
local class = require 'ext.class'
local SelfGrav = require 'solver.selfgrav'

-- TODO make this a behavior so Roe can be swapped with RoeImplicitLinear
local EulerRoe = class(SelfGrav(require 'solver.roe'))
--local EulerRoe = class(SelfGrav(require 'solver.roe_implicit_linearized'))

function EulerRoe:createEqn()
	self.eqn = require 'eqn.euler'(self)
end

return EulerRoe

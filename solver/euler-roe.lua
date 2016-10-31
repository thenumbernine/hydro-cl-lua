local class = require 'ext.class'
local Roe = require 'solver.roe'
local SelfGravitationBehavior = require 'solver.selfgrav'
local EulerRoe = class(SelfGravitationBehavior(Roe))
EulerRoe.eqn = require 'eqn/euler'()
return EulerRoe

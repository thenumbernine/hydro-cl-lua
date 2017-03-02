local class = require 'ext.class'
local Roe = require 'solver.roe'
local NoDiv = require 'solver.nodiv'
-- TODO - move NoDiv from solver to equation

local MHDRoe = class(NoDiv(Roe))

function MHDRoe:createEqn()
	self.eqn = require 'eqn.mhd'(self)
end

return MHDRoe

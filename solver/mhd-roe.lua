local class = require 'ext.class'
local Roe = require 'solver.roe'
local NoDiv = require 'solver.nodiv'

local MHDRoe = class(NoDiv(Roe))

function MHDRoe:createEqn()
	self.eqn = require 'eqn.mhd'(self)
end

return MHDRoe

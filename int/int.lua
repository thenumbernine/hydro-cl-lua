local class = require 'ext.class'
local real = require 'real'

local Integrator = class()

function Integrator:clearBuffer(buf)
	buf:fill(real(0))
end

return Integrator

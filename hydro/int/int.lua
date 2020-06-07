local class = require 'ext.class'
local real = require 'hydro.real'

local Integrator = class()

function Integrator:clearBuffer(buf)
	buf:fill(real(0))
end

return Integrator

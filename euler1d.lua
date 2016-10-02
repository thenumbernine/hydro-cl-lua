local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local Euler1D = class(Equation)

Euler1D.gamma = 7/5

Euler1D.consVars = table{'rho', 'mx', 'ETotal'}
Euler1D.primVars = table{'rho', 'vx', 'P'}
Euler1D.displayVars = table()
	:append(Euler1D.primVars)
	:append{'eInt', 'eKin', 'eTotal'} 

function Euler1D:header(clnumber)
	return '#define gamma '..clnumber(self.gamma)
end

function Euler1D:getTypeCode()
	return require 'makestruct'('prim_t', self.primVars) .. '\n'
		.. Euler1D.super.getTypeCode(self)
end

function Euler1D:solverCode()
	return '#include "euler1d.cl"'
end

-- TODO boundary methods, esp how to handle mirror

return Euler1D

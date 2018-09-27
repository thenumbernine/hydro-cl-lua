local table = require 'ext.table'
local class = require 'ext.class'
local Relaxation = require 'op.relaxation'

local Poisson = class(Relaxation)

Poisson.name = 'Poisson'

Poisson.solverCodeFile = 'op/poisson.cl'

function Poisson:getSolverCode()
	return table{
		Poisson.super.getSolverCode(self),
		self:getPoissonCode() or '',
	}:concat'\n'
end

return Poisson

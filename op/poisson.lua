local table = require 'ext.table'
local class = require 'ext.class'
local Relaxation = require 'op.relaxation'

local PoissonGaussSeidel = class(Relaxation)

PoissonGaussSeidel.name = 'PoissonGaussSeidel'

PoissonGaussSeidel.solverCodeFile = 'op/poisson.cl'

function PoissonGaussSeidel:getSolverCode()
	return table{
		PoissonGaussSeidel.super.getSolverCode(self),
		self:getPoissonCode() or '',
	}:concat'\n'
end

return PoissonGaussSeidel

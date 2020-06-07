--[[
from 2010 Colliander et al "Numerical Simulations ..."
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local GridSolver = require 'hydro.solver.gridsolver'

local NLSSolver = class(GridSolver)
NLSSolver.name = 'NonLinearSchrodinger'
NLSSolver.fixedDT = 1e-6
NLSSolver.eqnName = 'nls'

function NLSSolver:refreshCalcDTKernel() end
function NLSSolver:calcDT() return self.fixedDT end

return NLSSolver

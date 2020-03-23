--[[
based on 2014 Oliveira et al - Ergoregion instability- The hydrodynamic vortex
--]]
local class = require 'ext.class'
local GridSolver = require 'solver.gridsolver'

local WaveFDSolver = class(GridSolver)
WaveFDSolver.name = 'Wave_FiniteDifference'
WaveFDSolver.eqnName = 'wave-fd'
WaveFDSolver.fixedDT = 1e-5

-- TODO use the arbitrary finite difference stencil from bssn
WaveFDSolver.numGhost = 2

function WaveFDSolver:refreshCalcDTKernel() end
function WaveFDSolver:calcDT() 
	-- paper says h=1/30 ... 1/2000
	return .5 * tonumber(self.solverPtr.grid_dx.x)
	--return self.fixedDT 
end

return WaveFDSolver

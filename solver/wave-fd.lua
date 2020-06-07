--[[
based on 2014 Oliveira et al - Ergoregion instability- The hydrodynamic vortex
--]]
local class = require 'ext.class'
local GridSolver = require 'solver.gridsolver'

local WaveFDSolver = class(GridSolver)
WaveFDSolver.name = 'Wave-FD'
WaveFDSolver.eqnName = 'wave-fd'
WaveFDSolver.fixedDT = 1e-5

-- TODO use the arbitrary finite difference stencil from bssn
WaveFDSolver.numGhost = 2

function WaveFDSolver:init(...)
	WaveFDSolver.super.init(self, ...)
	-- fix the name
	self.name = WaveFDSolver.name..' '..self.integrator.name
end

function WaveFDSolver:refreshCalcDTKernel() end
function WaveFDSolver:calcDT() 
	-- paper says h=1/30 ... 1/2000
	return .5 * tonumber(self.solverPtr.grid_dx.x)
	--return self.fixedDT 
end


-- [[ here's all the work it takes to add 't' to a display var
-- how hard would it be to just put it in solver_t and update it immediately?
local real = require 'hydro.real'

WaveFDSolver.DisplayVar_U = class(WaveFDSolver.DisplayVar_U)

function WaveFDSolver.DisplayVar_U:setArgs(kernel)
	WaveFDSolver.DisplayVar_U.super.setArgs(self, kernel)
	kernel:setArg(4, real(self.solver.t))
end

function WaveFDSolver:getUBufDisplayVarsArgs()
	local args = WaveFDSolver.super.getUBufDisplayVarsArgs(self)
	args.extraArgs = args.extraArgs or {}
	table.insert(args.extraArgs, 'real t')
	args.extraArgNames = args.extraArgNames or {}
	table.insert(args.extraArgNames, 't')
	return args
end
--]]

return WaveFDSolver

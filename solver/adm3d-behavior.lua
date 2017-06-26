--[[
adds a constrainU call at the end of step()
seems the relativistic fluids have this as well
maybe I should move it to Solver?
--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local ADM3DEqn = require 'eqn.adm3d'

return function(parent)
	local template = class(parent)

	function template:createEqn()
		self.eqn = ADM3DEqn(self)
	end

	function template:refreshSolverProgram()
		template.super.refreshSolverProgram(self)

		self.constrainUKernel = self.solverProgram:kernel('constrainU', self.UBuf)
	end

	function template:step(dt)
		template.super.step(self, dt)

		self.app.cmds:enqueueNDRangeKernel{kernel=self.constrainUKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	end

	return template
end

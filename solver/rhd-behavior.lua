--[[
behavior for allocating a primbuf
and updating it using newton descent every iteration
and passing it as an extra arg to all those functions that need it

this is different from Euler in that Euler doesn't hold a prim buf
but maybe it should -- because this is running faster than Euler
--]]

local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'

return function(parent)
	local template = class(parent)

	--[[
	(TODO?) primBuf must be push/pop'd as well as UBuf
	(or am I safe just using the last iteration's prim values, and doing the newton descent to update them?)
	--]]
	function template:refreshSolverProgram()
		template.super.refreshSolverProgram(self)

		self.updatePrimsKernelObj = self.solverProgramObj:kernel('updatePrims', self.solverBuf, self.UBuf)
	end

	--[[ method 1: update prims after step() overall is called
	-- this leaves bad prims throughout RK4
	function template:step(dt)
		template.super.step(self, dt)

		self.updatePrimsKernelObj()
	end
	--]]
	-- [[ method 2: update the before calcDeriv
	-- this might take some extra calculations
	-- and prims could converge *from* different prevoius values even when converging *to* the same destination UBuf 
	function template:calcDeriv(derivBuf, dt)
		self.updatePrimsKernelObj()

		template.super.calcDeriv(self, derivBuf, dt)
	end
	--]]

	return template
end

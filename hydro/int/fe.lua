--[[
TODO 
instead of doing things like this
1) integrator clears deriv buf
2) solver writes to deriv buf
3) integrator adds deriv to U
it would be better to pass from integrator to solver
1) the final buffer to write
2) the write code

this could save us 1 less buffer, 1 less kernel call, and 1 less copy
--]]
local ffi = require 'ffi'
local CLBuffer = require 'cl.obj.buffer'
local roundup = require 'hydro.util.roundup'
local Integrator = require 'hydro.int.int'
local real = require 'hydro.real'

local ForwardEuler = Integrator:subclass()

ForwardEuler.name = 'forward Euler'

function ForwardEuler:init(solver)
	self.solver = solver
	self.derivBufObj = CLBuffer{
		env = solver.app.env,
		name = 'derivBuf',
		type = solver.eqn.symbols.cons_t,
		count = solver.numCells,
	}
end

function ForwardEuler:integrate(dt, callback)
	local solver = self.solver
	
	self:clearBuffer(self.derivBufObj)
	
	callback(self.derivBufObj)
	
	solver.multAddIntoKernelObj(solver.solverBuf, solver.UBuf, self.derivBufObj.obj, real(dt))

-- [[ I moved this from hydro/solver/gridsolver to integrator
-- this way I can do it after every substep in the RK integrator
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObj)) end
	solver:boundary()
if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
	solver:constrainU()
--]]
end

return ForwardEuler

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
local class = require 'ext.class'
local roundup = require 'util.roundup'
local Integrator = require 'int.int'
local CLBuffer = require 'cl.obj.buffer'

local ForwardEuler = class(Integrator)

ForwardEuler.name = 'forward Euler'

function ForwardEuler:init(solver)
	self.solver = solver
	self.derivBufObj = CLBuffer{
		env = solver.app.env,
		name = 'derivBuf',
		type = solver.eqn.cons_t,
		count = solver.numCells,
	}
end

local realptr = ffi.new'realparam[1]'
local function real(x)
	realptr[0] = x
	return realptr
end
function ForwardEuler:integrate(dt, callback)
	local solver = self.solver
	solver.app.cmds:enqueueFillBuffer{buffer=self.derivBufObj.obj, size=solver.numCells * ffi.sizeof(solver.eqn.cons_t)}
	callback(self.derivBufObj)
	solver.multAddKernelObj(solver.solverBuf, solver.UBuf, solver.UBuf, self.derivBufObj.obj, real(dt))
end

return ForwardEuler

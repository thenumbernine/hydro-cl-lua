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
local Integrator = require 'int.int'
local roundup = require 'roundup'

local ForwardEuler = class(Integrator)

ForwardEuler.name = 'forward Euler'

function ForwardEuler:init(solver)
	self.solver = solver
	self.derivBuf = solver.app.ctx:buffer{rw=true, size=solver.numCells * ffi.sizeof(solver.eqn.cons_t)}
end

local realptr = ffi.new'real[1]'
local function real(x)
	realptr[0] = x
	return realptr
end
function ForwardEuler:integrate(dt, callback)
	local solver = self.solver
	solver.app.cmds:enqueueFillBuffer{buffer=self.derivBuf, size=solver.numCells * ffi.sizeof(solver.eqn.cons_t)}
	callback(self.derivBuf)
	solver.multAddKernelObj(solver.UBuf, solver.UBuf, self.derivBuf, real(dt))
end

return ForwardEuler

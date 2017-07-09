local ffi = require 'ffi'
local class = require 'ext.class'
local Integrator = require 'int.int'
local roundup = require 'roundup'

local ForwardEuler = class(Integrator)

ForwardEuler.name = 'forward Euler'

function ForwardEuler:init(solver)
	self.solver = solver
	self.derivBuf = solver.app.ctx:buffer{rw=true, size=solver.volume * solver.eqn.numStates * ffi.sizeof(solver.app.real)}

	self.globalSize = roundup(solver.volume * solver.eqn.numStates, solver.maxWorkGroupSize)
end

function ForwardEuler:integrate(dt, callback)
	local solver = self.solver
	solver.app.cmds:enqueueFillBuffer{buffer=self.derivBuf, size=solver.volume * solver.eqn.numStates * ffi.sizeof(solver.app.real)}
	callback(self.derivBuf)
	solver.multAddKernel:setArgs(solver.UBuf, solver.UBuf, self.derivBuf, ffi.new('real[1]', dt))
	solver.app.cmds:enqueueNDRangeKernel{kernel=solver.multAddKernel, globalSize=self.globalSize, localSize=solver.localSize1d}
end

return ForwardEuler

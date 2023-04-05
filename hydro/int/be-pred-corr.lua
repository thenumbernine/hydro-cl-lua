--[[
this is gonna implement backward euler
under the assumption that I can turn it into the form:
https://en.wikipedia.org/wiki/Backward_Euler_method
the "non-stiff problems" implementation, which is assuming that the derivative calculation is not dependent on the state variable

U[0](t + Δt) = U(t)
U[i+1](t + Δt) = U(t) + Δt ∂U/∂t(t + Δt, U[i](t + Δt))
U(t + Δt) = lim i->∞ U[i](t + Δt)
--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local CLBuffer = require 'cl.obj.buffer'
local roundup = require 'hydro.util.roundup'
local Integrator = require 'hydro.int.int'
local real = require 'hydro.real'

local BackwardEulerPredCorr = class(Integrator)

BackwardEulerPredCorr.name = 'backward Euler predictor-corrector'

BackwardEulerPredCorr.maxiter = 10

function BackwardEulerPredCorr:init(solver, args)
	args = args or {}
	self.solver = solver
	self.maxiter = args.maxiter

	self.derivBufObj = CLBuffer{
		env = solver.app.env,
		name = 'derivBuf',
		type = solver.eqn.symbols.cons_t,
		count = solver.numCells,
	}

	self.U0BufObj = CLBuffer{
		env = solver.app.env,
		name = 'derivBuf',
		type = solver.eqn.symbols.cons_t,
		count = solver.numCells,
	}
end

function BackwardEulerPredCorr:integrate(dt, callback)
	local solver = self.solver
	local bufferSize = solver.numCells * ffi.sizeof(solver.eqn.symbols.cons_t)

	-- save U(t)
	solver.cmds:enqueueCopyBuffer{src=solver.UBuf, dst=self.U0BufObj.obj, size=bufferSize}
	
	for i=1,self.maxiter do
		-- calculate ∂U/∂t(t + Δt, U[i-1])
		self:clearBuffer(self.derivBufObj)
		callback(self.derivBufObj)
		
		-- U[i] = U(t) + Δt * ∂U/∂t(t + Δt, U[i-1])
		solver.multAddKernelObj(solver.solverBuf, solver.UBuf, self.U0BufObj.obj, self.derivBufObj.obj, real(dt))

		-- TODO move the checkNaNs into boundary() and constrainU() ?
		if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObj)) end
		solver:boundary()
		
		if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
		solver:constrainU()
	end
end

return BackwardEulerPredCorr

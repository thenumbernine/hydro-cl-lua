--[[
this is gonna implement backward euler
under the assumption that I can turn it into the form:
https://en.wikipedia.org/wiki/Backward_Euler_method
the "non-stiff problems" implementation, which is assuming that the derivative calculation is not dependent on the state variable

U[0](t + Δt) = U(t)
U[i+1](t + Δt) = U(t) + Δt ∂U/∂t(t + Δt, U[i](t + Δt))
U(t + Δt) = lim i→ ∞ U[i](t + Δt)
--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local math = require 'ext.math'
local CLBuffer = require 'cl.obj.buffer'
local roundup = require 'hydro.util.roundup'
local Integrator = require 'hydro.int.int'
local real = require 'hydro.real'

local half = require 'cl.obj.half'
local toreal, fromreal = half.toreal, half.fromreal

local BackwardEulerPredCorr = class(Integrator)

BackwardEulerPredCorr.name = 'backward Euler predictor-corrector'

BackwardEulerPredCorr.maxiter = 10

function BackwardEulerPredCorr:init(solver, args)
	args = args or {}
	self.solver = solver
	self.maxiter = args.maxiter or cmdline.intBEMaxIter
	self.epsilon = args.epsilon or cmdline.intBEEpsilon or nil
	self.verbose = cmdline.intVerbose or (args and args.verbose) or nil

	self.derivBufObj = CLBuffer{
		env = solver.app.env,
		name = 'derivBuf',
		type = solver.eqn.symbols.cons_t,
		count = solver.numCells,
	}

	-- iteration start state
	self.U0BufObj = CLBuffer{
		env = solver.app.env,
		name = 'U0Buf',
		type = solver.eqn.symbols.cons_t,
		count = solver.numCells,
	}

	if self.epsilon then
		-- previous iteration state
		self.U1BufObj = CLBuffer{
			env = solver.app.env,
			name = 'U1Buf',
			type = solver.eqn.symbols.cons_t,
			count = solver.numCells,
		}
	end
end

function BackwardEulerPredCorr:integrate(dt, callback)
	local solver = self.solver
	local bufferSize = solver.numCells * ffi.sizeof(solver.eqn.symbols.cons_t)

	-- save U(t)
	solver.cmds:enqueueCopyBuffer{src=solver.UBuf, dst=self.U0BufObj.obj, size=bufferSize}
	
	for i=1,self.maxiter do
		local prevUMag
		if self.epsilon then
			-- save U[i-1] for comparing later
			-- (todo? don't need ghost cells?)
			solver.cmds:enqueueCopyBuffer{src=solver.UBuf, dst=self.U1BufObj.obj, size=bufferSize}
			if cmdline.printBufsInt then
				print'\nU[i-1]:' solver:printBuf(solver.UBufObj, nil, solver.eqn.numIntStates)
			end
			-- calc its mag also
			solver.cmds:enqueueFillBuffer{buffer=solver.reduceBuf, size=solver.numCells * ffi.sizeof(solver.app.real)}
			solver.squareConsKernelObj.obj:setArg(2, solver.UBuf)
			solver.squareConsKernelObj()
			if cmdline.printBufsInt then
				print'\nreduceBuf <- |U[i-1]|^2:' solver:printBuf(solver.reduceBufObj)
			end
			prevUMag = math.sqrt(fromreal(solver.reduceSum(nil, solver.numCells)) / solver.volumeWithoutBorder)
			if self.verbose and type(self.verbose) == 'number' and self.verbose > 1 then
				print('prevUMag', prevUMag)
			end
		end

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
	
		-- optional epsilon
		-- subtract UBuf from previous UBuf and look at the magnitude
		-- compare to orig UBuf magnitude, break upon < epsilon
		if self.epsilon then
			solver.subtractIntoConsKernelObj(solver.solverBuf, self.U1BufObj.obj, solver.UBuf)
-- NOTICE this is clearing reduceBuf ... that's fine for calcDT and min across reduceBuf, right?
-- looks like in hydro/eqn/eqn.lua I"m writing dt=INFINITY even before reducing so yes.
-- better todo would be to just have dt not include ghost cells.
			solver.cmds:enqueueFillBuffer{buffer=solver.reduceBuf, size=solver.numCells * ffi.sizeof(solver.app.real)}
			solver.squareConsKernelObj.obj:setArg(2, self.U1BufObj.obj)
			solver.squareConsKernelObj()
			local deltaUMag = math.sqrt(fromreal(solver.reduceSum(nil, solver.numCells)) / solver.volumeWithoutBorder)
			if not math.isfinite(deltaUMag) then
				print("got non-finite residual: "..deltaUMag)
				return -- fail
			end
			if self.verbose and type(self.verbose) == 'number' and self.verbose > 1 then
				print('deltaUMag', deltaUMag)
			end
			local deltaRatio = deltaUMag/prevUMag
			-- TODO compare this to original
			if self.verbose then
				print('iter', i, 'res/prev umag: '..deltaRatio)
			end
			if deltaRatio < self.epsilon then break end
		end
	end
end

return BackwardEulerPredCorr

local ffi = require 'ffi'
local class = require 'ext.class'
local Integrator = require 'int.int'
local CLBuffer = require 'cl.obj.buffer'

local IterativeCrankNicolson = class(Integrator)
	
IterativeCrankNicolson.name = 'iterative Crank-Nicolson'

local realptr = ffi.new'realparam[1]'
local function real(x)
	realptr[0] = x
	return realptr
end	

IterativeCrankNicolson.iterations = 3

function IterativeCrankNicolson:init(solver)
	self.solver = solver
	self.iterations = args.iterations

	self.name = self.name..' '..self.iterations

	self.halfUBufObj = CLBuffer{
		env = solver.app.env,
		name = 'halfUBuf',
		type = solver.eqn.cons_t,
		count = solver.numCells,
	}

	self.derivBufObj = CLBuffer{
		env = solver.app.env,
		name = 'derivBuf',
		type = solver.eqn.cons_t,
		count = solver.numCells,
	}
end

-- u[i] = u[0] + dt * .5 * (du/dt(u[0]) + du/dt(u[i-1]))
function IterativeCrankNicolson:integrate(dt, callback)
	local solver = self.solver
	local bufferSize = solver.numCells * ffi.sizeof(solver.eqn.cons_t)
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[1])) end
		
	-- du/dt(u[0])
	callback(self.derivBufObj)

	-- u[.5] = u[0] + du/dt(u[0]) * .5 * dt
	solver.multAddKernelObj(solver.solverBuf, solver.halfUBuf, solver.UBuf, self.derivBufObj.obj, real(.5 * dt))

	-- u[1] = u[0] + du/dt(u[0]) * dt
	-- u[1] = u[.5] + du/dt(u[0]) * .5 * dt
	solver.multAddKernelObj(solver.solverBuf, solver.UBuf, solver.halfUBuf, self.derivBufObj.obj, real(.5 * dt))

	for i=1,self.iterations do
		-- du/dt(u[i-1])
		callback(self.derivBufObj)

		-- u[i] = u[0] + .5 * (du/dt(u[0]) + du/dt(u[i-1]))
		-- u[i] = u[.5] + du/dt(u[i-1]) * (.5 * dt)
		solver.multAddKernelObj(solver.solverBuf, solver.UBuf, solver.halfUBuf, self.derivBufObj.obj, real(.5 * dt))
	end
end

return IterativeCrankNicolson

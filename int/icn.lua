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

	self.srcUBufObj = CLBuffer{
		env = solver.app.env,
		name = 'srcUBuf',
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

--[[
k1 = .5 * f(t0, u0)
k2 = .5 * f(t0 + .5 * dt, u0 + .5 * dt * k1)
k3 = .5 * f(t0 + .5 * dt, u0 + .5 * dt * k2)
u_n+1 = u_n + dt * k3

I could implement this with a Butcher table using the 'RungeKutta' class,
but that would create 3 different buffers when I only need 1.
This makes me think, the RK class could be designed to reuse buffers that are no longer needed.
--]]

-- u[i] = u[0] + dt * .5 * (du/dt(u[0]) + du/dt(u[i-1]))
function IterativeCrankNicolson:integrate(dt, callback)
	local solver = self.solver
	local bufferSize = solver.numCells * ffi.sizeof(solver.eqn.cons_t)
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[1])) end
	
	-- TODO get rid of the cl.hpp impl in cl/*.lua and just replace it with cl/obj/*.lua
	solver.app.cmds:enqueueCopyBuffer{src=solver.UBuf, dst=self.srcUBufObj.obj, size=bufferSize}
	
	for i=1,self.iterations-1 do
		-- du/dt(u^i-1_k)
		callback(self.derivBufObj)

		-- u^i_k = u^0_k + du/dt(u^i-1_k) * .5 * dt
		solver.multAddKernelObj(
			solver.solverBuf,
			solver.UBuf,
			self.srcUBufObj.obj, 
			self.derivBufObj.obj,
			real(.5 * dt))
	end
		
	-- du/dt(u^i_k)
	callback(self.derivBufObj)

	-- u^m_k = u^0_k + du/dt(u^m-1_k) * dt
	solver.multAddKernelObj(
		solver.solverBuf,
		solver.UBuf,
		solver.srcUBufObj.obj,
		self.derivBufObj.obj,
		real(dt))
end

return IterativeCrankNicolson

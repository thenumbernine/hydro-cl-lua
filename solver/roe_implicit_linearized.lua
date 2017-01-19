local class = require 'ext.class'
local Roe = require 'solver.roe'
local ffi = require 'ffi'

local RoeImplicitLinearized = class(Roe)

function RoeImplicitLinearized:refreshGridSize(...)
	RoeImplicitLinearized.super.refreshGridSize(self, ...)

	-- the previous call in Solver:iterate is to self:calcDT
	-- which calls self.equation:calcInterfaceEigenBasis, which fills the self.eigenvalues table

	-- function that returns deriv when provided a state vector
	local function calc_dU_dt(dUdtBuf, UBuf)
		local oldUBuf = self.UBuf
		self.UBuf = UBuf
		self:calcDeriv(dUdtBuf, self.linearSolverDT)
		self.UBuf = oldUBuf
	end

	local linearSolverArgs = {
		--maxiter = 1000,
		x = self.UBuf,
		epsilon = 1e-10,
		maxiter = self.volume * self.eqn.numStates,
		restart = 10,
		-- logging:
		errorCallback = function(err, iter)
			print(self.t, iter, err)
		end,
	}

	-- [=[ backward Euler
	linearSolverArgs.b = self.lastUBuf
	linearSolverArgs.A = function(UNextBuf, UBuf)
		-- if this is a linearized implicit solver
		-- then the matrix should be computed before invoking the iterative solver
		-- which means the matrix coeffiicents shouldn't be changing per-iteration
		-- which means calc_dU_dt() should be based on the initial state and not the iterative state
		
		--UBuf = UBuf - dt * calc_dU_dt(self.UBuf)
		calc_dU_dt(dUdtbuf, self.UBuf)
		mulAdd(UNextBuf, UBuf, dUdtbuf, -self.linearSolverDT)
		
		-- ... but what about this?
		--UBuf = UBuf - dt * calc_dU_dt(UBuf)
		
		--self.boundaryMethod(UBuf)
	end
	--]=]
	--[=[ crank-nicolson - converges faster
	linearSolverArgs.b = (function(UBuf)
		UBuf = UBuf + .5 * dt * calc_dU_dt(UBuf)
		self.boundaryMethod(UBuf)
		return UBuf
	end)(UBuf)
	linearSolverArgs.A = function(UBuf)
		UBuf = UBuf - .5 * dt * calc_dU_dt(UBuf)
		self.boundaryMethod(UBuf)
		return UBuf
	end
	--]=]

	-- set up gmres solver here
	self.linearSolver = require 'solver.cl.gmres'(linearSolverArgs)
end

function RoeImplicitLinearized:createBuffers()
	RoeImplicitLinearized.super.createBuffers(self)
	self:clalloc('lastUBuf', self.volume * ffi.sizeof'cons_t')
end

-- step contains integrating flux and source terms
-- but not post iterate
function RoeImplicitLinearized:step(dt)
	self.linearSolverDT = dt
	self.linearSolver()
end

return RoeImplicitLinearized

local class = require 'ext.class'
local Roe = require 'solver.roe'

local RoeImplicitLinearized = class(Roe)

-- step contains integrating flux and source terms
-- but not post iterate
function RoeImplicitLinearized:step(dt)
	-- the previous call in Solver:iterate is to self:calcDT
	-- which calls self.equation:calcInterfaceEigenBasis, which fills the self.eigenvalues table

	-- function that returns deriv when provided a state vector
	local function calc_dU_dt(dUdtBuf, UBuf)
		local oldUBuf = self.UBuf
		self.UBuf = UBuf
		self:calcDeriv(dUdtBuf, dt)
		self.UBuf = oldUBuf
	end

-- [[ implicit via some linear solver
	local qs = self.qs

	local linearSolverArgs = {
		--maxiter = 1000,
		x = qs:clone(),
		epsilon = self.linearSolverEpsilon, 
		maxiter = self.linearSolverMaxIter,
		restart = self.linearSolverRestart,
		--[=[ mostly true ... mostly ...
		-- not true for any 2nd derivative terms
		-- this method is only used for Jacobi method, so I don't really care
		ADiag = (function()
			local n = self:newState()
			for i=1,self.gridsize do
				for j=1,self.numStates do
					n[i][j] = 1
				end
			end
			return n
		end)(),
		--]=]
		-- logging:
		errorCallback = self.errorLogging and function(err, convergenceIteration)
			print(self.t, convergenceIteration, err)
		end,
	}

	--[=[ identity.  do nothing.
	linearSolverArgs.b = qs:clone()
	linearSolverArgs.A = function(qs) return qs end
	--]=]
	-- [=[ backward Euler
	linearSolverArgs.b = qs:clone()
	linearSolverArgs.A = function(UNextBuf, UBuf)
		-- if this is a linearized implicit solver
		-- then the matrix should be computed before invoking the iterative solver
		-- which means the matrix coeffiicents shouldn't be changing per-iteration
		-- which means calc_dU_dt() should be based on the initial state and not the iterative state
		
		--UBuf = UBuf - dt * calc_dU_dt(self.UBuf)
		calc_dU_dt(dUdtbuf, self.UBuf)
		mulAdd(UNextBuf, UBuf, dUdtbuf, -dt)
		
		-- ... but what about this?
		--UBuf = UBuf - dt * calc_dU_dt(UBuf)
		
		--self.boundaryMethod(qs)
	end
	--]=]
	--[=[ crank-nicolson - converges faster
	linearSolverArgs.b = (function(qs)
		qs = qs + .5 * dt * calc_dU_dt(qs)
		self.boundaryMethod(qs)
		return qs
	end)(qs)
	linearSolverArgs.A = function(qs)
		qs = qs - .5 * dt * calc_dU_dt(qs)
		self.boundaryMethod(qs)
		return qs
	end
	--]=]

	self.qs = self.linearSolver(linearSolverArgs)
	if self.errorLogging then
		print()
	end
--]]
--[[ explicit - forward Euler - for debugging
	self.qs = self.qs + dt * calc_dU_dt(self.qs)
--]]
end

return RoeImplicitLinearized

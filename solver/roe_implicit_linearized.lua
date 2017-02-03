local class = require 'ext.class'
local table = require 'ext.table'
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

	-- TODO needs env
	-- for now I'll hack one into place from App and Solver
	-- but I should really rewrite this whole code to use cl.obj.env
	local fakeEnv = setmetatable({
		platform = self.app.platform,
		device = self.app.device,	
		ctx = self.app.ctx,
		cmds = self.app.cmds,
		totalGPUMem = 0,
		real = self.app.real,
		code = table{
			self.codePrefix,
			[[
//macro for the index
#define globalInt4()	(int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0)

//macros for arbitrary sizes
#define indexForInt4ForSize(i, sx, sy, sz) (i.x + sx * (i.y + sy * i.z))
#define initKernelForSize(sx, sy, sz) \
	int4 i = globalInt4(); \
	if (i.x >= sx || i.y >= sy || i.z >= sz) return; \
	int index = indexForInt4ForSize(i, sx, sy, sz);
]],
		}:concat'\n',
	}, require 'cl.obj.env')
	fakeEnv.base = fakeEnv:domain{
		size = require 'ext.range'(self.dim):map(function(i)
			return self.gridSize:ptr()[i+1]
		end),
	}
	local linearSolverArgs = {
		env = fakeEnv,
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
		
		-- use the integrator's deriv buffer and don't use this?
		local dUdtBuf = self.integrator.derivBuf
		
		--UBuf = UBuf - dt * calc_dU_dt(self.UBuf)
		calc_dU_dt(dUdtBuf, self.UBuf)
		self.linearSolver.args.mulAdd(UNextBuf, UBuf, dUdtBuf, -self.linearSolverDT)
		
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
	self:clalloc('lastUBuf', self.volume * ffi.sizeof(self.eqn.cons_t))
end

-- step contains integrating flux and source terms
-- but not post iterate
function RoeImplicitLinearized:step(dt)
	self.linearSolverDT = dt
	self.linearSolver()
end

return RoeImplicitLinearized

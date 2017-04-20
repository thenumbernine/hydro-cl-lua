local class = require 'ext.class'
local table = require 'ext.table'
local Roe = require 'solver.roe'
local ffi = require 'ffi'
local math = require 'ext.math'

local tolua = require 'ext.tolua'
local gputostr_n
local gputostr_cmds
function gputostr(x)
	local ptr = ffi.new('real['..gputostr_n..']')
	gputostr_cmds:enqueueReadBuffer{buffer=x, block=true, size=gputostr_n*ffi.sizeof'real', ptr=ptr}
	local s = ''
	for i=0,gputostr_n-1 do
		s = s .. '\t' ..ptr[i]
		if i % 6 == 5 then s = s .. '\n' end
	end
	return s 
end

local CLGMRES = require 'solver.cl.gmres'
local ThisGMRES = class(CLGMRES)

function ThisGMRES:newBuffer(name)
	if not self.cache then self.cache = {} end
	local cached = self.cache[name]
	if cached then return cached end
	cached = ThisGMRES.super.newBuffer(self, name)
	cached:fill()
	self.cache[name] = cached
	return cached
end

local RoeImplicitLinearized = class(Roe)

function RoeImplicitLinearized:refreshGridSize(...)
	RoeImplicitLinearized.super.refreshGridSize(self, ...)

	-- the previous call in Solver:iterate is to self:calcDT
	-- which calls self.equation:calcInterfaceEigenBasis, which fills the self.eigenvalues table

	-- function that returns deriv when provided a state vector
	-- dUdtBuf is a cl.buffer
	-- UBuf is a cl.obj.buffer
	local function calc_dU_dt(dUdtBuf, UBuf)
		-- setting self.UBuf won't make a difference because the kernels are already bound to the original self.UBuf
		-- copying over self.UBuf will overwrite the contents of 'x'
		-- so I need to (a1) make a new buffer for 'x'
		-- (a2) here, copy UBuf into self.UBuf before running the kernel
		-- or (b) have Solver:calcDeriv accept a param for the 'UBuf' 
		-- ... but that'll be hard since UBuf is exchanged with getULRBuf with usePLM
		self.UBufObj:copyFrom(UBuf)
		self.app.cmds:enqueueFillBuffer{buffer=dUdtBuf, size=self.volume * ffi.sizeof(self.eqn.cons_t)}
		self:calcDeriv(dUdtBuf, self.linearSolverDT)
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
gputostr_cmds = self.app.cmds
gputostr_n = self.volume * self.eqn.numStates

	local numreals = self.volume * self.eqn.numStates

	local function wrapBufferObj(field)
		return setmetatable({
			env = fakeEnv,
			name = field,
			type = 'real',
			size = self.volume * self.eqn.numStates,
			obj = self[field],
		}, require 'cl.obj.buffer')
	end

	self.UBufObj = wrapBufferObj'UBuf'
	self.RoeImplicitLinear_bObj = wrapBufferObj'RoeImplicitLinear_b'
	self.RoeImplicitLinear_xBufObj = wrapBufferObj'RoeImplicitLinear_x'

	local linearSolverArgs = {
		env = fakeEnv,
		--maxiter = 1000,
		x = self.RoeImplicitLinear_xBufObj,
		size = numreals,
		epsilon = 1e-10,
		maxiter = 10 * numreals,
		restart = 10,
		-- logging:
		errorCallback = function(err, iter, x, rLenSq, bLenSq)
--print('t', self.t, 'iter', iter, 'err', err, 'rLenSq', rLenSq, 'bLenSq', bLenSq)
			--print(' x', gputostr(x))
			if not math.isfinite(err) then
				error("got non-finite err: "..err)
			end
		end,
	}

	-- [=[ backward Euler
	linearSolverArgs.b = self.RoeImplicitLinear_bObj
	-- usage of A:
	-- A(r, x)
	-- A(w, v[i])
	linearSolverArgs.A = function(UNextBuf, UBuf)
		-- if this is a linearized implicit solver
		-- then the matrix should be computed before invoking the iterative solver
		-- which means the matrix coeffiicents shouldn't be changing per-iteration
		-- which means calc_dU_dt() should be based on the initial state and not the iterative state
--print('A begin UBuf', gputostr(UBuf.obj))
--print('A self.qs / RoeImplicitLinear_b', gputostr(self.RoeImplicitLinear_bObj.obj))
		-- use the integrator's deriv buffer and don't use this?
		local dUdtBuf = self.integrator.derivBuf
		calc_dU_dt(dUdtBuf, self.RoeImplicitLinear_bObj)
--print('A dq/dt / dUdtBuf', gputostr(dUdtBuf))
		
		--UNextBuf = UBuf - dt * calc_dU_dt(self.UBuf)
--print('A dt',self.linearSolverDT)		
		self.linearSolver.args.mulAdd(UNextBuf, UBuf, dUdtBuf, -self.linearSolverDT)
		
		--self.boundaryMethod(UBuf)
--print('A end UBuf', gputostr(UNextBuf.obj))		
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
	self.linearSolver = ThisGMRES(linearSolverArgs)
end

function RoeImplicitLinearized:createBuffers()
	RoeImplicitLinearized.super.createBuffers(self)
	self:clalloc('RoeImplicitLinear_b', self.volume * ffi.sizeof(self.eqn.cons_t))
	self:clalloc('RoeImplicitLinear_x', self.volume * ffi.sizeof(self.eqn.cons_t))
end

-- step contains integrating flux and source terms
-- but not post iterate
function RoeImplicitLinearized:step(dt)
--print('RoeImplicitLinearized:step() t',self.t,'dt',dt,'cfl',self.cfl,'gamma',self.eqn.guiVarsForName.heatCapacityRatio.value[0])
	self.linearSolverDT = dt
	
	-- UBuf needs to be overwritten to pass on to the calcFluxDeriv
	-- (TODO make calcFluxDeriv accept a parameter)
	self.RoeImplicitLinear_bObj:copyFrom(self.UBufObj)
	self.RoeImplicitLinear_xBufObj:copyFrom(self.UBufObj)
	self.linearSolver()
	self.UBufObj:copyFrom(self.RoeImplicitLinear_xBufObj)
--print('step done')
end

return RoeImplicitLinearized

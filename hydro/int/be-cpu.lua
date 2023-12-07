-- CPU version of the GMRES backward-Euler integrator

local ffi = require 'ffi'
local ig = require 'imgui'
local math = require 'ext.math'
local template = require 'template'
local CLBuffer = require 'cl.obj.buffer'
local Integrator = require 'hydro.int.int'

local BackwardEuler = Integrator:subclass()

BackwardEuler.name = 'backward Euler, CPU'

local matrix = require 'matrix'

function BackwardEuler:init(solver, args)
	self.solver = solver
	self.verbose = cmdline.intVerbose or (args and args.verbose) or nil

	-- gui vars:
	self.lastResidual = 0
	self.lastIter = 0

	-- full buffer, with ghost cells and non-integratable state variables
	self.derivBufObj = CLBuffer{
		env = solver.app.env,
		name = 'derivBuf',
		type = solver.app.real,
		count = solver.numCells * solver.eqn.numStates,
	}

	
	self.cpuUTemp = ffi.new(solver.eqn.symbols.cons_t..'[?]', solver.numCells)

	local volumeWithoutBorder = solver.volumeWithoutBorder
	local numRealsWithoutBorder = volumeWithoutBorder * solver.eqn.numIntStates
	self.krylov_b = matrix.zeros(numRealsWithoutBorder)
	self.krylov_zero = matrix.zeros(numRealsWithoutBorder)
end

function BackwardEuler:integrate(dt, callback)
	local solver = self.solver

	local function matToGpu(dstgpu, srcm)
		for k=0,tonumber(solver.sizeWithoutBorder.z-1) do
			for j=0,tonumber(solver.sizeWithoutBorder.y-1) do
				for i=0,tonumber(solver.sizeWithoutBorder.x-1) do
					for e=0,solver.eqn.numIntStates-1 do
						local dstIndex = 
							solver.stepSize.x * (solver.numGhost + i)
							+ solver.stepSize.y * ((solver.dim >= 2 and solver.numGhost or 0) + j)
							+ solver.stepSize.z * ((solver.dim >= 3 and solver.numGhost or 0) + k)
						local srcIndex = 1 + e + solver.eqn.numIntStates * tonumber(i + solver.sizeWithoutBorder.x * (j + solver.sizeWithoutBorder.y * k))
						self.cpuUTemp[dstIndex].ptr[e] = srcm[srcIndex]
					end
				end
			end
		end
		dstgpu:fromCPU(self.cpuUTemp)
	end

	local function gpuToMat(dstm, srcgpu)
		srcgpu:toCPU(self.cpuUTemp)
		for k=0,tonumber(solver.sizeWithoutBorder.z-1) do
			for j=0,tonumber(solver.sizeWithoutBorder.y-1) do
				for i=0,tonumber(solver.sizeWithoutBorder.x-1) do
					for e=0,solver.eqn.numIntStates-1 do
						local srcIndex = 
							solver.stepSize.x * (solver.numGhost + i)
							+ solver.stepSize.y * ((solver.dim >= 2 and solver.numGhost or 0) + j)
							+ solver.stepSize.z * ((solver.dim >= 3 and solver.numGhost or 0) + k)
						local dstIndex = 1 + e + solver.eqn.numIntStates * tonumber(i + solver.sizeWithoutBorder.x * (j + solver.sizeWithoutBorder.y * k))
						dstm[dstIndex] = self.cpuUTemp[srcIndex].ptr[e]
					end
				end
			end
		end
	end

	local volumeWithoutBorder = solver.volumeWithoutBorder
	local numRealsWithoutBorder = volumeWithoutBorder * solver.eqn.numIntStates

-- [[ evaluate derivative once, converge U to it
	self.derivBufObj:fill()
	callback(self.derivBufObj)
	local dU_dt = matrix(U)
	gpuToMat(dU_dt, self.derivBufObj)
--]]

	-- copy U(t) into 'b'
	gpuToMat(self.krylov_b, solver.UBufObj)
--print('krylov_b ... U(t)\n'..self.krylov_b)

	local restart = cmdline.intBERestart or (args and args.restart) or 20

	local UNext = require 'solver.gmres'{
		-- U(t+dt) = U(t) + dt dU/dt(t+dt)
		-- U(t+dt) = U(t) + dt F(U,t+dt)
		-- U(t+dt) - dt F(U,t+dt) = U(t)	<=> A(x) = b
		-- (I - dt dF/dU) U(t+dt) = U(t)
		-- so 'b' is U(t)
		-- and A(x) = x - dt * dU/dt(U,t+dt)
		A = function(U)
--[[ evaluate derivative while we converge state
			-- copy from U to solver.UBufObj (non-ghost-regions, only int vars)
			-- what about borders and non-int states?  
			--  ignore them and trust boundary() and constrainU() will deal with them
			matToGpu(solver.UBufObj, U)
			solver:boundary()
			solver:constrainU()
			self.derivBufObj:fill()
			callback(self.derivBufObj)
			
			local dU_dt = matrix(U)
			gpuToMat(dU_dt, self.derivBufObj)
			return self.krylov_b - dt * dU_dt
--]]
-- [[ evaluate derivative once, converte U to it
-- seems like this is just creating a linear transform that solves the forward-euler case ...
			return U - dt * dU_dt
--]]
--print('A(x) = U(t) - dt * dU/dt(t+dt)\n'..result)
		end,

		-- indistinguishable in the CPU version where the algorithm clones x initially
		b = self.krylov_b,
		
		errorCallback = function(residual, iter, x)
			self.lastResidual = residual
			self.lastIter = iter
			if self.verbose then
				print('t', solver.t, 'iter', iter, 'residual', residual)
			end
--if iter < numRealsWithoutBorder then return false end
			if not math.isfinite(residual) then
				print("got non-finite residual: "..residual)	-- error?
				return true	-- fail
			end
--			if residual < self.linearSolver.args.epsilon then return true end	
		end,
		
		clone = matrix,
	
		dot = function(a,b)
			return (a * b) / math.sqrt(numRealsWithoutBorder)
		end,
		
		zero = self.krylov_zero,
		epsilon = cmdline.intBEEpsilon or (args and args.epsilon) or 1e-10,
		restart = restart,
		maxiter = cmdline.intBEMaxIter or (args and args.maxiter) or restart * numRealsWithoutBorder,
	}

	matToGpu(solver.UBufObj, UNext)
end

return BackwardEuler

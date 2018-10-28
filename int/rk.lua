local ffi = require 'ffi'
local class = require 'ext.class'
local Integrator = require 'int.int'

local RungeKutta = class(Integrator)

function RungeKutta:init(solver)
	self.solver = solver
	self.order = #self.alphas
	assert(#self.betas == self.order)
	for i=1,self.order do
		assert(#self.alphas[i] == self.order)
		assert(#self.betas[i] == self.order)
	end

	self.UBufs = {}
	self.derivBufs = {}
	for i=1,self.order do
		local needed = false
		for m=i,self.order do
			needed = needed or self.alphas[m][i] ~= 0
		end
		if needed then
			self.UBufs[i] = solver.app.ctx:buffer{rw=true, size=solver.numCells * ffi.sizeof(solver.eqn.cons_t)}
		end
	
		local needed = false
		for m=i,self.order do
			needed = needed or self.betas[m][i] ~= 0
		end
		if needed then
			self.derivBufs[i] = solver.app.ctx:buffer{rw=true, size=solver.numCells * ffi.sizeof(solver.eqn.cons_t)}
		end
	end
end

local realptr = ffi.new'realparam[1]'
local function real(x)
	realptr[0] = x
	return realptr
end
function RungeKutta:integrate(dt, callback)
	local solver = self.solver
	local bufferSize = solver.numCells * ffi.sizeof(solver.eqn.cons_t)
	
	solver.multAddKernelObj.obj:setArgs(solver.solverBuf, solver.UBuf, solver.UBuf)
	
	--u(0) = u^n
	local needed = false
	for m=1,self.order do
		needed = needed or self.alphas[m][1] ~= 0
	end
	if needed then
		solver.app.cmds:enqueueCopyBuffer{src=solver.UBuf, dst=self.UBufs[1], size=bufferSize}
	end

	--L(u^(0))
	local needed = false
	for m=1,self.order do
		needed = needed or self.betas[m][1] ~= 0
	end
	if needed then
		solver.app.cmds:enqueueFillBuffer{buffer=self.derivBufs[1], size=bufferSize}
	end

	for i=2,self.order+1 do
		--u^(i) = sum k=0 to i-1 of (alpha_ik u^(k) + dt beta_ik L(u^(k)) )
		solver.app.cmds:enqueueFillBuffer{buffer=solver.UBuf, size=bufferSize}
		for k=1,i-1 do
			if self.alphas[i-1][k] ~= 0 then
				solver.multAddKernelObj.obj:setArg(3, self.UBufs[k])
				solver.multAddKernelObj.obj:setArg(4, real(self.alphas[i-1][k]))
				solver.multAddKernelObj()
			end
			if self.betas[i-1][k] ~= 0 then
				solver.multAddKernelObj.obj:setArg(3, self.derivBufs[k])
				solver.multAddKernelObj.obj:setArg(4, real(self.betas[i-1][k] * dt))
				solver.multAddKernelObj()
			end
		end
	
		if i <= self.order then
			--only do this if alpha_mi != 0 for any m
			--otherwise there's no need to store this buffer
			local needed = false
			for m=i,self.order do
				needed = needed or self.alphas[m][i] ~= 0
			end
			if needed then
				solver.app.cmds:enqueueCopyBuffer{src=solver.UBuf, dst=self.UBufs[i], size=bufferSize}
			end
		
			--likewise here, only if beta_mi != 0 for any m
			--with that in mind, no need to allocate these buffers unless they are needed.
			local needed = false
			for m=i,self.order do
				needed = needed or self.betas[m][i] ~= 0
			end
			if needed then
				solver.app.cmds:enqueueFillBuffer{buffer=self.derivBufs[i], size=bufferSize}
				callback(self.derivBufs[i])
			end
		end
		--else just leave the state in there
	end
end

return RungeKutta

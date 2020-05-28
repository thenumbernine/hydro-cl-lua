local ffi = require 'ffi'
local class = require 'ext.class'
local Integrator = require 'int.int'
local CLBuffer = require 'cl.obj.buffer'
local real = require 'real'

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
	self.derivBufObjs = {}
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
			self.derivBufObjs[i] = CLBuffer{
				env = solver.app.env,
				name = 'derivBuf'..i,
				type = solver.eqn.cons_t,
				count = solver.numCells,
			}
			
			self:clearBuffer(self.derivBufObjs[i])

if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[i])) end
		end
	end
end

function RungeKutta:integrate(dt, callback)
	local solver = self.solver
	local bufferSize = solver.numCells * ffi.sizeof(solver.eqn.cons_t)
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[1])) end

--print("integrating rk with dt "..dt)
--print('order = '..self.order)	
	solver.multAddKernelObj.obj:setArgs(solver.solverBuf, solver.UBuf, solver.UBuf)
	
	--u(0) = u^n
	local needed = false
	for m=1,self.order do
		needed = needed or self.alphas[m][1] ~= 0
	end
	if needed then
--print('UBufs[1] = UBuf')
		solver.cmds:enqueueCopyBuffer{src=solver.UBuf, dst=self.UBufs[1], size=bufferSize}
	end

	--L(u^(0))
	local needed = false
	for m=1,self.order do
		needed = needed or self.betas[m][1] ~= 0
	end
	if needed then
--print('derivBufObjs[1].obj = dU/dt(UBuf)')
if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[1])) end
	
		self:clearBuffer(self.derivBufObjs[1])

if solver.checkNaNs then solver.cmds:finish() assert(solver:checkFinite(self.derivBufObjs[1])) end
		callback(self.derivBufObjs[1])
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[1])) end
	
if cmdline.printBufs then
	print()
	print('UBuf RK4 deriv 1:')
	solver:printBuf(self.derivBufObjs[1])
end

	-- else error"Why do you have a RK butcher table without using the first iteration derivative?"
	end

	for i=2,self.order+1 do
		--u^(i) = sum k=0 to i-1 of (alpha_ik u^(k) + dt beta_ik L(u^(k)) )
--io.write('UBuf = 0')	
		self:clearBuffer(solver.UBufObj)
		for k=1,i-1 do
--io.write(' + UBufs['..k..'] * (a['..(i-1)..']['..k..'] = '..self.alphas[i-1][k]..')')
			if self.alphas[i-1][k] ~= 0 then
				solver.multAddKernelObj.obj:setArg(3, self.UBufs[k])
				solver.multAddKernelObj.obj:setArg(4, real(self.alphas[i-1][k]))
				solver.multAddKernelObj()
if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
			end
		end
		for k=1,i-1 do
--io.write(' + derivBufObjs['..k..'] * (b['..(i-1)..']['..k..'] = '..self.betas[i-1][k]..') * dt')
			if self.betas[i-1][k] ~= 0 then
				solver.multAddKernelObj.obj:setArg(3, self.derivBufObjs[k].obj)
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[k])) end
				solver.multAddKernelObj.obj:setArg(4, real(self.betas[i-1][k] * dt))
				solver.multAddKernelObj()
if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
			end
		end
--print()

-- [[ I moved this from solver/gridsolver to integrator
-- this way I can do it after every substep in the RK integrator
if solver.checkNaNs then assert(solver:checkFinite(derivBufObj)) end
		solver:boundary()
if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
		solver:constrainU()
--]]

if cmdline.printBufs then
	print()
	print('UBuf after RK4 step'..(i-1)..':')
	solver:printBuf(solver.UBufObj)
end
		if i <= self.order then
			--only do this if alpha_mi != 0 for any m
			--otherwise there's no need to store this buffer
			local needed = false
			for m=i,self.order do
				needed = needed or self.alphas[m][i] ~= 0
			end
			if needed then
--print('UBufs['..i..'] = UBuf')
				solver.cmds:enqueueCopyBuffer{src=solver.UBuf, dst=self.UBufs[i], size=bufferSize}
			end
		
			--likewise here, only if beta_mi != 0 for any m
			--with that in mind, no need to allocate these buffers unless they are needed.
			local needed = false
			for m=i,self.order do
				needed = needed or self.betas[m][i] ~= 0
			end
			if needed then
--print('derivBufObjs['..i..'].obj = dU/dt(UBuf)')				
				
				self:clearBuffer(self.derivBufObjs[i])

if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[i])) end
				callback(self.derivBufObjs[i])
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[i])) end

if cmdline.printBufs then
	print()
	print('UBuf RK4 deriv '..i..':')
	solver:printBuf(self.derivBufObjs[i])
end
			end
		end
		--else just leave the state in there	
	end

end

return RungeKutta

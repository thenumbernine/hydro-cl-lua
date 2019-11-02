local ffi = require 'ffi'
local class = require 'ext.class'
local Integrator = require 'int.int'
local CLBuffer = require 'cl.obj.buffer'

local RungeKutta = class(Integrator)

local realptr = ffi.new'realparam[1]'
local function real(x)
	realptr[0] = x
	return realptr
end
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
			
			self.derivBufObjs[i]:fill(real(0))

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
		solver.app.cmds:enqueueCopyBuffer{src=solver.UBuf, dst=self.UBufs[1], size=bufferSize}
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
	
		self.derivBufObjs[1]:fill(real(0))

if solver.checkNaNs then solver.app.cmds:finish() assert(solver:checkFinite(self.derivBufObjs[1])) end
		callback(self.derivBufObjs[1])
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[1])) end
	
if cmdline.printBufs then
	print('UBuf RK4 deriv 1:')
	solver:printBuf(self.derivBufObjs[1], nil, nil, 24)
end

	-- else error"Why do you have a RK butcher table without using the first iteration derivative?"
	end

	for i=2,self.order+1 do
		--u^(i) = sum k=0 to i-1 of (alpha_ik u^(k) + dt beta_ik L(u^(k)) )
--io.write('UBuf = 0')	
		solver.UBufObj:fill(real(0))
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
		if solver.useConstrainU then
			solver.constrainUKernelObj()
if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
		end
--]]

if cmdline.printBufs then
	print('UBuf after RK4 step'..(i-1)..':')
	solver:printBuf(solver.UBufObj, nil, nil, 25)
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
				solver.app.cmds:enqueueCopyBuffer{src=solver.UBuf, dst=self.UBufs[i], size=bufferSize}
			end
		
			--likewise here, only if beta_mi != 0 for any m
			--with that in mind, no need to allocate these buffers unless they are needed.
			local needed = false
			for m=i,self.order do
				needed = needed or self.betas[m][i] ~= 0
			end
			if needed then
--print('derivBufObjs['..i..'].obj = dU/dt(UBuf)')				
				
				self.derivBufObjs[i]:fill(real(0))

if solver.checkNaNs then assert(solver:checkFinite(solver.UBufObj)) end
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[i])) end
				callback(self.derivBufObjs[i])
if solver.checkNaNs then assert(solver:checkFinite(self.derivBufObjs[i])) end

if cmdline.printBufs then
	print('UBuf RK4 deriv '..i..':')
	solver:printBuf(self.derivBufObjs[i], nil, nil, 24)
end
			end
		end
		--else just leave the state in there	
	end

end

return RungeKutta

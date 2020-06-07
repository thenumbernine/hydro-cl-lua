#!/usr/bin/env luajit
require 'ext'
local ffi = require 'ffi'
local unistd = require 'ffi.c.unistd'
require 'ffi.c.stdlib'
local dirp = unistd.getcwd(nil, 0)
local dir = ffi.string(dirp)
ffi.C.free(dirp)
unistd.chdir'../..'

-- honestly what's App used for anyways, beyond the gui?
-- the cl.env ... can I build a solver without app, but just with a cl.env?
local App = class(require 'hydro.app')

function App:setup()
	local args = {
		app = self, 
		dim = 1,
		integrator = 'forward Euler',
		fluxLimiter = 'superbee',
		coord = 'cartesian',
		mins = {-1, -1, -1},
		maxs = {1, 1, 1},
		gridSize = {256, 256, 256},
		boundary = {
			xmin='freeflow',
			xmax='freeflow',
			ymin='freeflow',
			ymax='freeflow',
			zmin='freeflow',
			zmax='freeflow',
		},
		initState = 'Sod',
		--initState = 'self-gravitation test 1',
	}

	-- [=[ two-solver testing ...
	-- running two solvers at once causes errors
	local app = self
	local cl = require 'hydro.solver.roe' 
	args.eqn = 'euler'
	--args.eqn = 'maxwell'
	--args.eqn = 'mhd'
	self.solvers:insert(cl(args))
	self.solvers:insert(cl(args))
	
	local s1, s2 = self.solvers:unpack()

	local numReals = s1.numCells * s1.eqn.numStates
	local ptr1 = ffi.new('real[?]', numReals)
	local ptr2 = ffi.new('real[?]', numReals)
	local function compare(buf1, buf2)
		app.cmds:enqueueReadBuffer{buffer=buf1, block=true, size=ffi.sizeof(app.real) * numReals, ptr=ptr1}
		app.cmds:enqueueReadBuffer{buffer=buf2, block=true, size=ffi.sizeof(app.real) * numReals, ptr=ptr2}
		local diff
		for i=0,numReals-1 do
			if ptr1[i] ~= ptr2[i] then
				diff = i 
				break
			end
		end
		if diff then
			for i=0,numReals-1 do
				io.write(('%7d '):format(i))
				for j,v in ipairs{ptr1, ptr2} do
					io.write(({' ','/'})[j])
					local vi = v[i]
					local s = vi and ('%.8e'):format(v[i]) or 'nil'
					local ip = ffi.cast('int*', v+i)
					s = s .. '(' .. (ip and (('%08x'):format(ip[1]):sub(-8) .. ('%08x'):format(ip[0]):sub(-8)) or 'nil') .. ')'
					io.write(s)
				end
				local col = 1
				if i%col==col-1 then print() end
			end
			print()
			local ch = diff % s1.eqn.numStates
			local x = math.floor(diff / tonumber(s1.eqn.numStates * s1.gridSize.x)) % tonumber(s1.gridSize.y)
			local y = math.floor(diff / tonumber(s1.eqn.numStates * s1.gridSize.x * s1.gridSize.y))
			print('index '..diff
				..' coord '..x..', '..y..' ch '..ch
				..' differs:',ptr1[diff], ptr2[diff])
			s1:save's1'
			s2:save's2'
			error'here'
		end
	end

	--for _,s in ipairs(self.solvers) do s.useFixedDT = true s.fixedDT = .025 end
	--s2.integrator = s1.integrator
	--function s2:update() self.t = self.t + .1 end
	--function s2:boundary() end
	--function s2:calcDT() return s1.fixedDT end
	-- [==[ messing with step to find the problem
	function s2:step(dt)
		--[[ fails
		s1.integrator:integrate(dt, function(derivBufObj)
			self:calcDeriv(derivBufObj, dt)
		end)
		--]]
		--[[ seems to work fine, but we're not adding to UBuf
		self:calcDeriv(s1.integrator.derivBufObj, dt)
		--]]
		--[[ works as well .. but doesn't add to s2's UBuf
		self:calcDeriv(self.integrator.derivBufObj, dt)
		--]]
		-- [[ fails with inline forward euler
		-- calcDeriv runs fine on its own
		-- everything except calcDeriv runs on its own
		-- but as soon as the two are put together, it dies
		app.cmds:enqueueFillBuffer{buffer=self.integrator.derivBufObj.obj, size=self.numCells * self.eqn.numStates * ffi.sizeof(app.real)}
		-- this produces crap in s2
		-- if it's not added into s2's UBuf then we're safe
		-- if it isn't called and zero is added to s2's UBuf then we're safe
		self:calcDeriv(self.integrator.derivBufObj, dt)
		
	-- why does derivBufObj differ?
	-- something is corrupting every numStates data ... like something is writing out of bounds ...
	--compare(s1.integrator.derivBufObj.obj, s2.integrator.derivBufObj.obj)
		
		self.multAddKernelObj(
			self.UBuf,
			s1.UBuf,
			
			-- using self.integrator.derivBufObj.obj fails
			-- but using s1.integrator.derivBufObj.obj works ... worked ....
			s1.integrator.derivBufObj.obj,
			
			ffi.new('real[1]', dt))
		--]]
		--[[ fails with the other solver's forward-euler and multAddKernel
		app.cmds:finish()
		app.cmds:enqueueFillBuffer{buffer=s1.integrator.derivBufObj.obj, size=self.numCells * self.eqn.numStates * ffi.sizeof(app.real)}
		self:calcDeriv(s1.integrator.derivBufObj, dt)
		s1.multAddKernelObj(self.UBuf, self.UBuf, s1.integrator.derivBufObj.obj, ffi.new('real[1]', dt))
		app.cmds:finish()
		--]]
		--[[ just adding s1's deriv to s2?  works fine
		s1.multAddKernelObj(self.UBuf, self.UBuf, s1.integrator.derivBufObj.obj, ffi.new('real[1]', dt))
		s1.multAddKernelObj.obj:setArgs(s1.UBuf, s1.UBuf, s1.integrator.derivBufObj.obj, ffi.new('real[1]', dt))
		--]]
		-- so the code in common is when calcDeriv is called by the 2nd solver ...
		-- ... regardless of what buffer it is written to
	end
	--]==]
	--[==[ comparing buffers.  tends to die on the boundaries even if it is working (why is that?)
	function s2:update()
		s1.update(self)
		-- ...annd even when using s1's derivBufObj, this dies once the wave hits a boundary
		-- complains about negative'd values (with mirror boundary conditions)
		compare(s1.UBuf, s2.UBuf)
	end
	--]==]
	--]=]
end

App():run()

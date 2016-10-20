local ffi = require 'ffi'
local gl = require 'ffi.OpenGL'
local class = require 'ext.class'
local string = require 'ext.string'
local table = require 'ext.table'
local range = require 'ext.range'
local vec3sz = require 'vec3sz'
local vec3 = require 'vec.vec3'
local clnumber = require 'clnumber'



local Integrator = class()


local ForwardEuler = class(Integrator)
ForwardEuler.name = 'forward Euler'

function ForwardEuler:init(solver)
	self.solver = solver
	self.derivBuf = solver.app.ctx:buffer{rw=true, size=solver.volume * solver.eqn.numStates * ffi.sizeof(solver.app.real)}
end

function ForwardEuler:integrate(dt, callback)
	local solver = self.solver
	callback(self.derivBuf)
	solver.multAddKernel:setArgs(solver.UBuf, solver.UBuf, self.derivBuf, ffi.new('real[1]', dt))
	solver.app.cmds:enqueueNDRangeKernel{kernel=solver.multAddKernel, globalSize=solver.volume * solver.eqn.numStates, localSize=solver.localSize1d}
end

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
			self.UBufs[i] = solver.app.ctx:buffer{rw=true, size=solver.volume * solver.eqn.numStates * ffi.sizeof(solver.app.real)}
		end
	
		local needed = false
		for m=i,self.order do
			needed = needed or self.betas[m][i] ~= 0
		end
		if needed then
			self.derivBufs[i] = solver.app.ctx:buffer{rw=true, size=solver.volume * solver.eqn.numStates * ffi.sizeof(solver.app.real)}
		end
	end
end

function RungeKutta:integrate(dt, callback)
	local solver = self.solver
	local realSize = ffi.sizeof(solver.app.real)
	local length = solver.volume * solver.eqn.numStates
	local bufferSize = length * realSize
	
	solver.multAddKernel:setArgs(solver.UBuf, solver.UBuf)
	
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
				solver.multAddKernel:setArg(2, self.UBufs[k])
				solver.multAddKernel:setArg(3, ffi.new('real[1]', self.alphas[i-1][k]))
				solver.app.cmds:enqueueNDRangeKernel{kernel=solver.multAddKernel, globalSize=length, localSize=solver.localSize1d}
			end
			if self.betas[i-1][k] ~= 0 then
				solver.multAddKernel:setArg(2, self.derivBufs[k])
				solver.multAddKernel:setArg(3, ffi.new('real[1]', self.betas[i-1][k] * dt))
				solver.app.cmds:enqueueNDRangeKernel{kernel=solver.multAddKernel, globalSize=length, localSize=solver.localSize1d}
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

--the following are from https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Classic_fourth-order_method

local RungeKutta2 = class(RungeKutta)
RungeKutta2.name = 'Runge-Kutta 2'
RungeKutta2.alphas = {
	{1,0},
	{1,0},
}
RungeKutta2.betas = {
	{.5, 0},
	{0, 1},
}

local RungeKutta2Heun = class(RungeKutta)
RungeKutta2Heun.name = 'Runge-Kutta 2 Heun'
RungeKutta2Heun.alphas = { 
	{1, 0},
	{1, 0},
}
RungeKutta2Heun.betas = {
	{1, 0},
	{.5, .5},
}

local RungeKutta2Ralston = class(RungeKutta)
RungeKutta2Ralston.name = 'Runge-Kutta 2 Ralston'
RungeKutta2Ralston.alphas = {
	{1, 0},
	{1, 0},
}
RungeKutta2Ralston.betas = {
	{2./3., 0},
	{1./4., 3./4.},
}

local RungeKutta3 = class(RungeKutta)
RungeKutta3.name = 'Runge-Kutta 3'
RungeKutta3.alphas = {
	{1, 0, 0},
	{1, 0, 0},
	{1, 0, 0},
}
RungeKutta3.betas = {
	{.5, 0, 0},
	{-1, 2, 0},
	{1./6., 2./6., 1./6.},
}

local RungeKutta4 = class(RungeKutta)
RungeKutta4.name = 'Runge-Kutta 4'
RungeKutta4.alphas = {
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{1, 0, 0, 0},
}
RungeKutta4.betas = {
	{.5, 0, 0, 0},
	{0, .5, 0, 0},
	{0, 0, 1, 0},
	{1./6., 2./6., 2./6., 1./6.},
}

local RungeKutta4_3_8thsRule = class(RungeKutta)
RungeKutta4_3_8thsRule.name = 'Runge-Kutta 4, 3/8ths rule'
RungeKutta4_3_8thsRule.alphas = {
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{1, 0, 0, 0},
}
RungeKutta4_3_8thsRule.betas = {
	{1./3., 0, 0, 0},
	{-1./3., 0, 0, 0},
	{1, -1, 1, 0},
	{1./8., 3./8., 3./8., 1./8.},
}

--the following are from http://www.ams.org/journals/mcom/1998-67-221/S0025-5718-98-00913-2/S0025-5718-98-00913-2.pdf

local RungeKutta2TVD = class(RungeKutta)
RungeKutta2TVD.name = 'Runge-Kutta 2, TVD'
RungeKutta2TVD.alphas = {
	{1, 0},
	{.5, .5},
}
RungeKutta2TVD.betas = {
	{1, 0},
	{0, .5},
}

local RungeKutta2NonTVD = class(RungeKutta)
RungeKutta2NonTVD.name = 'Runge-Kutta 2, non-TVD'
RungeKutta2NonTVD.alphas = {
	{1, 0},
	{1, 0},
}
RungeKutta2NonTVD.betas = {
	{-20, 0},
	{41./40., -1./40.},
}

local RungeKutta3TVD = class(RungeKutta)
RungeKutta3TVD.name = 'Runge-Kutta 3, TVD'
RungeKutta3TVD.alphas = {
	{1, 0, 0},
	{3./4., 1./4., 0},
	{1./3., 2./3., 0},
}
RungeKutta3TVD.betas = {
	{1, 0, 0},
	{0, 1./4., 0},
	{0, 0, 2./3.},
}

local RungeKutta4TVD = class(RungeKutta)
RungeKutta4TVD.name = 'Runge-Kutta 4, TVD'
RungeKutta4TVD.alphas = {
	{1, 0, 0, 0},
	{649./1600., 951./1600., 0, 0},
	{53989./2500000., 4806213./20000000., 23619./32000., 0},
	{1./5., 6127./30000., 7873./30000., 1./3.},
}
RungeKutta4TVD.betas = {
	{.5, 0, 0, 0},
	{-10890423./25193600., 5000./7873, 0, 0},
	{-102261./5000000., -5121./20000., 7873./10000., 0},
	{1./10., 1./6., 0, 1./6.},
}

-- this one is from http://lsec.cc.ac.cn/lcfd/DEWENO/paper/WENO_1996.pdf

local RungeKutta4NonTVD = class(RungeKutta)
RungeKutta4NonTVD.name = 'Runge-Kutta 4, non-TVD'
RungeKutta4NonTVD.alphas = {
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{1, 0, 0, 0},
	{-1./3., 1./3., 2./3., 1./3.},
}
RungeKutta4NonTVD.betas = {
	{.5, 0, 0, 0},
	{0, .5, 0, 0},
	{0, 0, 1, 0},
	{0, 0, 0, 1./6.},
}

local xs = table{'x', 'y', 'z'}
local minmaxs = table{'min', 'max'}

local Solver = class()

Solver.name = 'Roe'
Solver.numGhost = 2

Solver.integrators = {
	ForwardEuler,
	RungeKutta2,
	RungeKutta2Heun,
	RungeKutta2Ralston,
	RungeKutta3,
	RungeKutta4,
	RungeKutta4_3_8thsRule,
	RungeKutta2TVD,
	RungeKutta2NonTVD,
	RungeKutta3TVD,
	RungeKutta4TVD,
	RungeKutta4NonTVD,
}
Solver.integratorNames = table.map(Solver.integrators, function(integrator) return integrator.name end)

--[[
args:
	app
	gridSize
	dim (optional, required if gridSize is a vec3sz)
--]]
function Solver:init(args)
	assert(args)
	local gridSize = assert(args.gridSize)
	if type(gridSize) == 'number' then 
		self.gridSize = vec3sz(gridSize,1,1)
		assert(not args.dim or args.dim == 1)
		self.dim = 1
	elseif type(gridSize) == 'cdata' 
	and ffi.typeof(gridSize) == vec3sz 
	then
		self.gridSize = vec3sz(gridSize)
		self.dim = assert(args.dim)
	elseif type(gridSize) == 'table' then
		self.gridSize = vec3sz(table.unpack(gridSize))
		assert(not args.dim or args.dim <= #gridSize)
		self.dim = args.dim or #gridSize
	else
		error("can't understand args.gridSize type "..type(args.gridSize).." value "..tostring(args.gridSize))
	end

	for i=self.dim,2 do self.gridSize:ptr()[i] = 1 end

	self.color = vec3(math.random(), math.random(), math.random()):normalize()

	self.mins = vec3(table.unpack(args.mins or {-1, -1, -1}))
	self.maxs = vec3(table.unpack(args.maxs or {1, 1, 1}))

	self.t = 0
	
	self.eqn = assert(args.eqn)
	self.name = self.eqn.name..' '..self.name

	self.app = assert(args.app)
	
	-- https://stackoverflow.com/questions/15912668/ideal-global-local-work-group-sizes-opencl
	-- product of all local sizes must be <= max workgroup size
	local maxWorkGroupSize = tonumber(self.app.device:getInfo'CL_DEVICE_MAX_WORK_GROUP_SIZE')

	self.offset = vec3sz(0,0,0)
	self.localSize1d = maxWorkGroupSize 
	
--	self.localSize = self.dim < 3 and vec3sz(16,16,16) or vec3sz(4,4,4)
print('maxWorkGroupSize',maxWorkGroupSize)
	self.localSize = vec3sz(1,1,1)
	local rest = maxWorkGroupSize
	local localSizeX = math.min(gridSize[1], 2^math.ceil(math.log(rest^(1/self.dim),2)))
	self.localSize.x = localSizeX
	rest = rest / localSizeX
	if self.dim == 2 then
		self.localSize.y = rest
		self.localSize.z = 1
	elseif self.dim == 3 then
		local localSizeY = math.min(gridSize[2], 2^math.ceil(math.log(math.sqrt(rest),2)))
		self.localSize.y = localSizeY
		local localSizeZ = rest / localSizeY
		self.localSize.z = localSizeZ
	end
print('self.localSize',self.localSize)

	self.useFixedDT = ffi.new('bool[1]', false)
	self.fixedDT = ffi.new('float[1]', 0)
	self.cfl = ffi.new('float[1]', .5)
	
	self.initStatePtr = ffi.new('int[1]', (table.find(self.eqn.initStateNames, args.initState) or 1)-1)
	
	self.integratorPtr = ffi.new('int[1]', (self.integratorNames:find(args.integrator) or 1)-1)
	
	self.slopeLimiterPtr = ffi.new('int[1]', (self.app.slopeLimiterNames:find(args.slopeLimiter) or 1)-1)

	self.boundaryMethods = {}
	for i=1,3 do
		for _,minmax in ipairs(minmaxs) do
			local var = xs[i]..minmax
			self.boundaryMethods[var] = ffi.new('int[1]', self.app.boundaryMethods:find(
					(args.boundary or {})[var] or 'freeflow'
				)-1)
		end
	end

	self:refreshGridSize()
end

function Solver:refreshGridSize()

	self.volume = tonumber(self.gridSize:volume())
	self.dxs = vec3(range(3):map(function(i)
		return (self.maxs[i] - self.mins[i]) / tonumber(self.gridSize:ptr()[i-1])	
	end):unpack())

	self:refreshIntegrator()	-- depends on eqn & gridSize

	self:createDisplayVars()	-- depends on eqn
	self:createBuffers()		-- depends on eqn & gridSize
	self:createCodePrefix()		-- depends on eqn, gridSize, displayVars

	self:refreshInitStateProgram()
	self:refreshSolverProgram()
	self:refreshBoundaryProgram()

	self:resetState()
end

function Solver:refreshIntegrator()
	self.integrator = self.integrators[self.integratorPtr[0]+1](self)
end

function Solver:createDisplayVars()
	self.displayVars = table()
	local function makevars(buffer, ...)
		for i=1,select('#',...) do
			local name = tostring(select(i,...))
			self.displayVars:insert{
				buffer = buffer,
				name = buffer..'_'..name,
				enabled = ffi.new('bool[1]', 
					buffer == 'U' and (self.dim==1 or i==1)
					or (buffer == 'error' and self.dim==1)
				),
				color = vec3(math.random(), math.random(), math.random()):normalize(),
				--heatMapTexPtr = ffi.new('int[1]', 0),	-- hsv, isobar, etc ...
				heatMapFixedRangePtr = ffi.new('bool[1]', true),
				heatMapValueMinPtr = ffi.new('float[1]', 0),
				heatMapValueMaxPtr = ffi.new('float[1]', 1),
			}
		end
	end
	
	makevars('U', table.unpack(self.eqn.displayVars))
	makevars('wave', range(0,self.eqn.numWaves-1):unpack())
	makevars('eigen', table.unpack(self.eqn:getEigenInfo().displayVars))
	makevars('deltaUTilde', range(0,self.eqn.numWaves-1):unpack())
	makevars('rTilde', range(0,self.eqn.numWaves-1):unpack())
	makevars('flux', range(0,self.eqn.numStates-1):unpack())
	makevars('reduce', '0')	-- might contain nonsense :-p
	makevars('deriv', range(0,self.eqn.numStates-1):unpack())
	makevars('error', 'ortho', 'flux')
end

local errorType = 'error_t'
local errorTypeCode = 'typedef struct { real ortho, flux; } '..errorType..';'

function Solver:createBuffers()
	local ctx = self.app.ctx
	local realSize = ffi.sizeof(self.app.real)

	-- to get sizeof
	ffi.cdef(self.eqn:getEigenInfo().typeCode)
	
	ffi.cdef(errorTypeCode)

	self.UBuf = ctx:buffer{rw=true, size=self.volume * self.eqn.numStates * realSize}
	self.reduceBuf = ctx:buffer{rw=true, size=self.volume * realSize}
	self.reduceResultPtr = ffi.new('real[1]', 0)
	self.reduceSwapBuf = ctx:buffer{rw=true, size=self.volume * realSize / self.localSize1d}
	self.waveBuf = ctx:buffer{rw=true, size=self.volume * self.dim * self.eqn.numWaves * realSize}
	self.eigenBuf = ctx:buffer{rw=true, size=self.volume * self.dim * ffi.sizeof'eigen_t'}
	self.deltaUTildeBuf = ctx:buffer{rw=true, size=self.volume * self.dim * self.eqn.numWaves * realSize}
	self.rTildeBuf = ctx:buffer{rw=true, size=self.volume * self.dim * self.eqn.numWaves * realSize}
	self.fluxBuf = ctx:buffer{rw=true, size=self.volume * self.dim * self.eqn.numStates * realSize}
	
	-- debug only
	self.fluxXformBuf = ctx:buffer{rw=true, size=self.volume * self.dim * ffi.sizeof'fluxXform_t'}
	local errorTypeSize = ffi.sizeof(errorType)
	self.errorBuf = ctx:buffer{rw=true, size=self.volume * self.dim * errorTypeSize}

	-- CL/GL interop

	local GLTex2D = require 'gl.tex2d'
	self.tex = GLTex2D{
		width = self.gridSize.x,
		height = self.gridSize.y,
		internalFormat = gl.GL_RGBA32F,
		format = gl.GL_RGBA,
		type = gl.GL_FLOAT,
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_LINEAR,
		wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT},
	}

	if self.app.useGLSharing then
		local ImageGL = require 'cl.imagegl'
		self.texCLMem = ImageGL{context=ctx, tex=self.tex, write=true}
	else
		self.calcDisplayVarToTexPtr = ffi.new(self.app.real..'[?]', self.volume)
		
		--[[ PBOs?
		self.calcDisplayVarToTexPBO = ffi.new('gl_int[1]', 0)
		gl.glGenBuffers(1, self.calcDisplayVarToTexPBO)
		gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, self.calcDisplayVarToTexPBO[0])
		gl.glBufferData(gl.GL_PIXEL_UNPACK_BUFFER, self.tex.width * self.tex.height * ffi.sizeof(self.app.real) * 4, nil, gl.GL_STREAM_READ)
		gl.glBindBuffer(gl.GL_PIXEL_UNPACK_BUFFER, 0)
		--]]
	end
end

function Solver:createCodePrefix()
	local lines = table()
	if self.dim == 3 then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable'
	end
	if self.app.real == 'double' then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_fp64 : enable'
	end

	lines:append(table{'',2,4,8}:map(function(n)
		return 'typedef '..self.app.real..n..' real'..n..';'
	end)):append{
		'#define dim '..self.dim,
		'#define numGhost '..self.numGhost,
		'#define numStates '..self.eqn.numStates,
		'#define numWaves '..self.eqn.numWaves,
	}:append(xs:map(function(x,i)
		return '#define mins_'..x..' '..clnumber(self.mins[i])..'\n'
			.. '#define maxs_'..x..' '..clnumber(self.maxs[i])..'\n'
	end)):append{
		'constant real4 mins = (real4)(mins_x, '..(self.dim<2 and '0' or 'mins_y')..', '..(self.dim<3 and '0' or 'mins_z')..', 0);', 
		'constant real4 maxs = (real4)(maxs_x, '..(self.dim<2 and '0' or 'maxs_y')..', '..(self.dim<3 and '0' or 'maxs_z')..', 0);', 
	}:append(xs:map(function(x,i)
		return '#define d'..x..' '..clnumber(self.dxs[i])
	end)):append{
		'constant real4 dxs = (real4)(dx, dy, dz, 0);',
		'#define dx_min '..clnumber(math.min(table.unpack(self.dxs, 1, self.dim))),
	}:append(xs:map(function(name,i)
		return '#define gridSize_'..name..' '..tonumber(self.gridSize[name])
	end)):append{
		'constant int4 gridSize = (int4)(gridSize_x, gridSize_y, gridSize_z, 0);',
		'constant int4 stepsize = (int4)(1, gridSize_x, gridSize_x * gridSize_y, gridSize_x * gridSize_y * gridSize_z);',
		'#define INDEX(a,b,c)	((a) + gridSize_x * ((b) + gridSize_y * (c)))',
		'#define INDEXV(i)		INDEX((i).x, (i).y, (i).z)',
		'#define CELL_X(i) (real4)('
			..'(real)(i.x + .5) * dx + mins_x, '
			..(--self.dim < 2 and '0,' or 
				'(real)(i.y + .5) * dy + mins_y, '
			)
			..(--self.dim < 3 and '0,' or 
				'(real)(i.z + .5) * dz + mins_z, '
			)
			..'0);',
	}:append{
		self.eqn.getTypeCode and self.eqn:getTypeCode() or nil
	}:append{
		-- run here for the code, and in buffer for the sizeof()
		self.eqn:getEigenInfo().typeCode
	}:append{
		-- bounds-check macro
		'#define OOB(lhs,rhs) (i.x < lhs || i.x >= gridSize_x - rhs'
			.. (self.dim < 2 and '' or ' || i.y < lhs || i.y >= gridSize_y - rhs')
			.. (self.dim < 3 and '' or ' || i.z < lhs || i.z >= gridSize_z - rhs')
			.. ')',
		-- define i, index, and bounds-check
		'#define SETBOUNDS(lhs,rhs)	\\',
		'\tint4 i = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0); \\',
		'\tif (OOB(lhs,rhs)) return; \\',
		'\tint index = INDEXV(i);',
	}

	lines:append(self.displayVars:map(function(var,i)
		return '#define display_'..var.name..' '..i
	end))

	-- output the first and last indexes of display vars associated with each buffer
	local buffers = table()
	for _,var in ipairs(self.displayVars) do
		buffers[var.buffer] = buffers[var.buffer] or table()
		buffers[var.buffer]:insert(var)
	end
	for buffer,vars in pairs(buffers) do
		lines:insert('#define displayFirst_'..buffer..' display_'..vars[1].name)
		lines:insert('#define displayLast_'..buffer..' display_'..vars:last().name)
	end

	self.codePrefix = lines:concat'\n'
end

function Solver:refreshInitStateProgram()
	local initStateCode = table{
		self.codePrefix,
		self.eqn:getInitStateCode(self, clnumber),
	}:concat'\n'
	self.initStateProgram = require 'cl.program'{context=self.app.ctx, devices={self.app.device}, code=initStateCode}
	self.initStateKernel = self.initStateProgram:kernel('initState', self.UBuf)
end

function Solver:resetState()
	self.app.cmds:finish()
	self.app.cmds:enqueueNDRangeKernel{kernel=self.initStateKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.app.cmds:finish()
end

-- depends on buffers
function Solver:refreshSolverProgram()

	-- code that depend on real and nothing else
	-- TODO move to app, along with reduceBuf


	local commonCode = table():append
		{self.app.is64bit and '#pragma OPENCL EXTENSION cl_khr_fp64 : enable' or nil
	}:append{
		'typedef '..self.app.real..' real;',
	
		--templates in C ...
		'#define reduce_accum_init INFINITY',
		'#define reduce_operation(x,y) min(x,y)',
		'#define reduce_name reduceMin',
		'#include "reduce.cl"',
		'#undef reduce_accum_init',
		'#undef reduce_operation',
		'#undef reduce_name',

		'#define reduce_accum_init -INFINITY',
		'#define reduce_operation(x,y) max(x,y)',
		'#define reduce_name reduceMax',
		'#include "reduce.cl"',
		'#undef reduce_accum_init',
		'#undef reduce_operation',
		'#undef reduce_name',
	
		[[
__kernel void multAdd(
	__global real* a,
	const __global real* b,
	const __global real* c,
	real d
) {
	size_t i = get_global_id(0);
	if (i >= get_global_size(0)) return;
	a[i] = b[i] + c[i] * d;
}	
]],
	}:concat'\n'

	self.commonProgram = require 'cl.program'{context=self.app.ctx, devices={self.app.device}, code=commonCode}

	for _,name in ipairs{'Min', 'Max'} do
		self['reduce'..name..'Kernel'] = self.commonProgram:kernel(
			'reduce'..name,
			self.reduceBuf,
			{ptr=nil, size=self.localSize1d * ffi.sizeof(self.app.real)},
			ffi.new('int[1]', self.volume),
			self.reduceSwapBuf)
	end
	
	self.multAddKernel = self.commonProgram:kernel'multAdd'


	-- solver code


	local slopeLimiterCode = 'real slopeLimiter(real r) {'
		.. self.app.slopeLimiters[1+self.slopeLimiterPtr[0]].code 
		.. '}'
		
	local code = table{
		self.codePrefix,
		slopeLimiterCode,
		self.eqn:getEigenInfo().code or '',
		
		'typedef struct { real min, max; } range_t;',
		errorTypeCode,	
		self.eqn:solverCode(self) or '',
	}:append(self.app.useGLSharing and {
		'#define calcDisplayVar_dstImage_t '..(self.dim == 3 and 'image3d_t' or 'image2d_t'),
		'#define calcDisplayVar_writeImageArgs '..(self.dim == 3 and '(int4)(i.x, i.y, i.z, 0)' or '(int2)(i.x, i.y)'),
		
		'#define calcDisplayVar_name calcDisplayVarToTex',
		'#define calcDisplayVar_output_tex',
		'#include "calcDisplayVar.cl"',
		'#undef calcDisplayVar_name',
		'#undef calcDisplayVar_output_tex',
	} or {}):append{	
		'#define calcDisplayVar_name calcDisplayVarToBuffer',
		'#define calcDisplayVar_output_buffer',
		'#include "calcDisplayVar.cl"',
		'#undef calcDisplayVar_name',
		'#undef calcDisplayVar_output_buffer',

		'#include "solver.cl"',
	}:concat'\n'

	self.solverProgram = require 'cl.program'{context=self.app.ctx, code=code}
	local success, message = self.solverProgram:build{self.app.device}
	if not success then
		-- show code
		print(string.split(string.trim(code),'\n'):map(function(l,i) return i..':'..l end):concat'\n')
		-- show errors
--		message = string.split(string.trim(message),'\n'):filter(function(line) return line:find'error' end):concat'\n'
		error(message)	
	end

	self.calcDTKernel = self.solverProgram:kernel('calcDT', self.reduceBuf, self.UBuf);
	
	self.calcEigenBasisKernel = self.solverProgram:kernel('calcEigenBasis', self.waveBuf, self.eigenBuf, self.fluxXformBuf, self.UBuf)
	self.calcDeltaUTildeKernel = self.solverProgram:kernel('calcDeltaUTilde', self.deltaUTildeBuf, self.UBuf, self.eigenBuf)
	self.calcRTildeKernel = self.solverProgram:kernel('calcRTilde', self.rTildeBuf, self.deltaUTildeBuf, self.waveBuf)
	self.calcFluxKernel = self.solverProgram:kernel('calcFlux', self.fluxBuf, self.UBuf, self.waveBuf, self.eigenBuf, self.deltaUTildeBuf, self.rTildeBuf)
	
	self.calcDerivFromFluxKernel = self.solverProgram:kernel'calcDerivFromFlux'
	self.calcDerivFromFluxKernel:setArg(1, self.fluxBuf)
	if self.eqn.useSourceTerm then
		self.addSourceTermKernel = self.solverProgram:kernel'addSourceTerm'
		self.addSourceTermKernel:setArg(1, self.UBuf)
	end

	if self.app.useGLSharing then
		self.calcDisplayVarToTexKernel = self.solverProgram:kernel('calcDisplayVarToTex', self.texCLMem)
	end
	self.calcDisplayVarToBufferKernel = self.solverProgram:kernel('calcDisplayVarToBuffer', self.reduceBuf)

	self.calcErrorsKernel = self.solverProgram:kernel('calcErrors', self.errorBuf, self.waveBuf, self.eigenBuf, self.fluxXformBuf)	
end

function Solver:refreshBoundaryProgram()
	local lines = table()
	lines:insert(self.codePrefix)
	lines:insert[[
__kernel void boundary(
	__global cons_t* UBuf
) {
]]
	if self.dim == 1 then 
		lines:insert[[
	if (get_global_id(0) != 0) return;
]]
	elseif self.dim == 2 then
		lines:insert[[
	int i = get_global_id(0);
]]
	elseif self.dim == 3 then
		lines:insert[[
	int2 i = (int2)(get_global_id(0), get_global_id(1));
]]
	end
	
	-- 1D: use a small 1D kernel and just run once 
	-- 2D: use a 1D kernel the size of the max dim 
		
	lines:insert[[
	for (int j = 0; j < numGhost; ++j) {
]]

	for side=1,self.dim do

		if self.dim == 2 then
			if side == 1 then
				lines:insert[[
	if (i < gridSize_y) {
]]
			elseif side == 2 then
				lines:insert[[
	if (i < gridSize_x) {
]]
			end
		elseif self.dim == 3 then
			if side == 1 then
				lines:insert[[
	if (i.x < gridSize_y && i.y < gridSize_z) {
]]
			elseif side == 2 then
				lines:insert[[
	if (i.x < gridSize_x && i.y < gridSize_z) {
]]
			elseif side == 3 then
				lines:insert[[
	if (i.x < gridSize_x && i.y < gridSize_y) {
]]
			end
		end

		local function index(j)
			if self.dim == 1 then
				return j
			elseif self.dim == 2 then
				if side == 1 then
					return 'INDEX('..j..',i,0)'
				elseif side == 2 then
					return 'INDEX(i,'..j..',0)'
				end
			elseif self.dim == 3 then
				if side == 1 then
					return 'INDEX('..j..',i.x,i.y)'
				elseif side == 2 then
					return 'INDEX(i.x,'..j..',i.y)'
				elseif side == 3 then
					return 'INDEX(i.x,i.y,'..j..')'
				end
			end
			error'TODO'
		end

		local x = xs[side]

		lines:insert(({
			periodic = '\t\tUBuf['..index'j'..'] = UBuf['..index'gridSize_x-2*numGhost+j'..'];',
			mirror = table{
				'\t\tUBuf['..index'j'..'] = UBuf['..index'2*numGhost-1-j'..'];',
			}:append(table.map((self.eqn.mirrorVars or {})[side] or {}, function(var)
				return '\t\tUBuf['..index'j'..'].'..var..' = -UBuf['..index'j'..'].'..var..';'
			end)):concat'\n',
			freeflow = '\t\tUBuf['..index'j'..'] = UBuf['..index'numGhost'..'];',
		})[self.app.boundaryMethods[1+self.boundaryMethods[x..'min'][0]]])

		lines:insert(({
			periodic = '\t\tUBuf['..index'gridSize_x-numGhost+j'..'] = UBuf['..index'numGhost+j'..'];',
			mirror = table{
				'\t\tUBuf['..index'gridSize_x-numGhost+j'..'] = UBuf['..index'gridSize_x-numGhost-1-j'..'];',
			}:append(table.map((self.eqn.mirrorVars or {})[side] or {}, function(var)
				return '\t\tUBuf['..index'gridSize_x-numGhost+j'..'].'..var..' = -UBuf['..index'gridSize_x-numGhost+j'..'].'..var..';'
			end)):concat'\n',
			freeflow = '\t\tUBuf['..index'gridSize_x-numGhost+j'..'] = UBuf['..index'gridSize_x-numGhost-1'..'];',
		})[self.app.boundaryMethods[1+self.boundaryMethods[x..'max'][0]]])
	
		if self.dim == 2 then
			lines:insert'\t}'
		end
	end
	
	lines:insert'\t}'
	lines:insert'}'

	local code = lines:concat'\n'
	
	self.app.cmds:finish()
	
	self.boundaryProgram = require 'cl.program'{context=self.app.ctx, devices={self.app.device}, code=code} 
	self.boundaryKernel = self.boundaryProgram:kernel('boundary', self.UBuf);
end

function Solver:boundary()
	-- 1D:
	if self.dim == 1 then
		self.app.cmds:enqueueNDRangeKernel{kernel=self.boundaryKernel, globalSize=self.localSize1d, localSize=self.localSize1d}
	elseif self.dim == 2 then
		local maxSize = math.max(tonumber(self.gridSize.x), tonumber(self.gridSize.y))
		self.app.cmds:enqueueNDRangeKernel{kernel=self.boundaryKernel, globalSize=maxSize, localSize=self.localSize1d}
	elseif self.dim == 3 then
	else
		error("can't run boundary for dim "..tonumber(self.dim))
	end
end

function Solver:reduceMin()
	return self:reduce(self.reduceMinKernel)
end

function Solver:reduceMax()
	return self:reduce(self.reduceMaxKernel)
end

function Solver:reduce(kernel)
	local reduceSize = self.volume
	local dst = self.reduceSwapBuf
	local src = self.reduceBuf
	local iter = 0

	local reduceLocalSize1D = 16
	
	while reduceSize > 1 do
		iter = iter + 1
		--TODO instead of >> 4, make sure it matches whatever reduceLocalSize1D is
		-- ... which just so happens to be 16 (i.e. 1 << 4) at the moment
		local nextSize = bit.rshift(reduceSize, 4)
		if 0 ~= bit.band(reduceSize, bit.lshift(1, 4) - 1) then 
			nextSize = nextSize + 1 
		end
		local reduceGlobalSize = math.max(reduceSize, reduceLocalSize1D)
		kernel:setArg(0, src)
		kernel:setArg(2, ffi.new('int[1]', reduceSize))
		kernel:setArg(3, dst)
		self.app.cmds:enqueueNDRangeKernel{kernel=kernel, dim=1, globalSize=reduceGlobalSize, localSize=math.min(reduceGlobalSize, reduceLocalSize1D)}
		--self.app.cmds:finish()
		dst, src = src, dst
		reduceSize = nextSize
	end
	self.app.cmds:enqueueReadBuffer{buffer=src, block=true, size=ffi.sizeof(self.app.real), ptr=self.reduceResultPtr}
	return self.reduceResultPtr[0]
end

function Solver:calcDeriv(derivBuf, dt)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcEigenBasisKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcErrorsKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDeltaUTildeKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcRTildeKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.calcFluxKernel:setArg(6, ffi.new('real[1]', dt))
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcFluxKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}

	-- calcDerivFromFlux zeroes the derivative buffer
	self.calcDerivFromFluxKernel:setArg(0, derivBuf)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDerivFromFluxKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}

	-- addSourceTerm adds to the derivative buffer
	if self.eqn.useSourceTerm then
		self.addSourceTermKernel:setArg(0, derivBuf)
		self.app.cmds:enqueueNDRangeKernel{kernel=self.addSourceTermKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	end
end

function Solver:integrate(derivBuf, dt)
	-- forward Euler
	self.multAddKernel:setArgs(self.UBuf, self.UBuf, derivBuf, ffi.new('real[1]', dt))
	self.app.cmds:enqueueNDRangeKernel{kernel=self.multAddKernel, globalSize=self.volume * self.eqn.numStates, localSize=self.localSize1d}
end

function Solver:update()
	self:boundary()

	-- calc cell wavespeeds -> dts
	if self.useFixedDT[0] then
		dt = tonumber(self.fixedDT[0])
	else
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDTKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		dt = tonumber(self.cfl[0]) * self:reduceMin()
	end
	self.fixedDT[0] = dt
	
	self.integrator:integrate(dt, function(derivBuf)
		self:calcDeriv(derivBuf, dt)
	end)

	self.t = self.t + dt
end

function Solver:calcDisplayVarToTex(varIndex)
	local var = self.displayVars[varIndex]
	if self.app.useGLSharing then
		-- copy to GL using cl_*_gl_sharing
		gl.glFinish()
		self.app.cmds:enqueueAcquireGLObjects{objs={self.texCLMem}}
		self.calcDisplayVarToTexKernel:setArg(1, ffi.new('int[1]', varIndex))
		self.calcDisplayVarToTexKernel:setArg(2, self[var.buffer..'Buf'])
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDisplayVarToTexKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		self.app.cmds:enqueueReleaseGLObjects{objs={self.texCLMem}}
		self.app.cmds:finish()
	else
		-- download to CPU then upload with glTexSubImage2D
		local ptr = self.calcDisplayVarToTexPtr
		local tex = self.tex
		self.calcDisplayVarToBufferKernel:setArg(1, ffi.new('int[1]', varIndex))
		self.calcDisplayVarToBufferKernel:setArg(2, self[var.buffer..'Buf'])
		self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDisplayVarToBufferKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
		self.app.cmds:enqueueReadBuffer{buffer=self.reduceBuf, block=true, size=ffi.sizeof(self.app.real) * self.volume, ptr=ptr}
		if self.app.is64bit then
			for i=0,self.volume-1 do
				ffi.cast('float*',ptr)[i] = ffi.cast('double*',ptr)[i]
			end
		end
		tex:bind()
		gl.glTexSubImage2D(gl.GL_TEXTURE_2D, 0, 0, 0, tex.width, tex.height, gl.GL_RED, gl.GL_FLOAT, ptr)
		tex:unbind()
		require 'gl.report' 'here'
	end
end

function Solver:calcDisplayVarRange(varIndex)
	local var = self.displayVars[varIndex]
	self.calcDisplayVarToBufferKernel:setArg(1, ffi.new('int[1]', varIndex))
	self.calcDisplayVarToBufferKernel:setArg(2, self[var.buffer..'Buf'])
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDisplayVarToBufferKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	local ymin = self:reduceMin()
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDisplayVarToBufferKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	local ymax = self:reduceMax()
	return ymin, ymax
end

return Solver

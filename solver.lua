local ffi = require 'ffi'
local gl = require 'ffi.OpenGL'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local vec3sz = require 'vec3sz'
local vec3 = require 'vec.vec3'

local xs = table{'x', 'y', 'z'}
local minmaxs = table{'min', 'max'}

local Solver = class()

Solver.name = 'Roe'
Solver.numGhost = 2

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
		error("can't understand args.gridSize type "..type(args.gridsize).." value "..tostring(args.gridSize))
	end

	for i=self.dim,2 do self.gridSize:ptr()[i] = 1 end

	self.color = vec3(math.random(), math.random(), math.random()):normalize()
	self.volume = tonumber(self.gridSize:volume())

	self.mins = args.mins or {-1, -1, -1}
	self.maxs = args.maxs or {1, 1, 1}
	self.t = 0
	
	self.eqn = assert(args.eqn)
	self.name = self.eqn.name..' '..self.name

	self.offset = vec3sz(0,0,0)
	self.localSize1d = 16
	self.localSize = self.dim < 3 and vec3sz(16,16,16) or vec3sz(8,8,8)

	self.app = assert(args.app)

	self.useFixedDT = ffi.new('bool[1]', false)
	self.fixedDT = ffi.new('float[1]', 0)
	self.cfl = ffi.new('float[1]', .5)
	
	self.initState = ffi.new('int[1]', (table.find(self.eqn.initStates, args.initState) or 1)-1)
	self.slopeLimiter = ffi.new('int[1]', (self.app.slopeLimiterNames:find(args.slopeLimiter) or 1)-1)

	self.boundaryMethods = {}
	for i=1,self.dim do
		for _,minmax in ipairs(minmaxs) do
			local var = xs[i]..minmax
			self.boundaryMethods[var] = ffi.new('int[1]', self.app.boundaryMethods:find'freeflow'-1)
		end
	end

	self:refreshGridSize()
end

local function clnumber(x)
	local s = tostring(x)
	if s:find'e' then return s end
	if not s:find('%.') then s = s .. '.' end
	return s
end

function Solver:refreshGridSize()

	self:createDisplayVars()	-- depends on eqn
	self:createBuffers()		-- depends on eqn & gridsize
	self:createCodePrefix()		-- depends on eqn, gridsize, displayvars

	self:refreshSolverProgram()
	self:refreshBoundaryProgram()

	self:resetState()
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
					buffer == 'U'
					or buffer == 'fluxError'
					or buffer == 'orthoError'
				),
				color = vec3(math.random(), math.random(), math.random()):normalize(),
			}
		end
	end
	
	makevars('U', self.eqn.displayVars:unpack())
	makevars('wave', range(0,self.eqn.numWaves-1):unpack())
	
	-- TODO should I allow this?
	-- it restricts the eigen struct to only contain reals ...
--	makevars('eigen', range(0,self.eqn.numEigen-1):unpack())
	
	makevars('dt', '0')
	makevars('deltaUTilde', range(0,self.eqn.numWaves-1):unpack())
	makevars('rTilde', range(0,self.eqn.numWaves-1):unpack())
	makevars('flux', range(0,self.eqn.numStates-1):unpack())
	makevars('deriv', range(0,self.eqn.numStates-1):unpack())
	makevars('orthoError', '0')
	makevars('fluxError', '0')
end

function Solver:createBuffers()
	local ctx = self.app.ctx
	local realSize = ffi.sizeof(self.app.real)

	ffi.cdef(self.eqn:getEigenTypeCode())

	self.UBuf = ctx:buffer{rw=true, size=self.volume * self.eqn.numStates * realSize}
	self.reduceBuf = ctx:buffer{rw=true, size=self.volume * realSize}
	self.reduceResultPtr = ffi.new('real[1]', 0)
	self.reduceSwapBuf = ctx:buffer{rw=true, size=self.volume * realSize / self.localSize1d}
	self.waveBuf = ctx:buffer{rw=true, size=self.volume * self.dim * self.eqn.numWaves * realSize}
	self.eigenBuf = ctx:buffer{rw=true, size=self.volume * self.dim * ffi.sizeof(self.eqn.eigenType)}
	self.deltaUTildeBuf = ctx:buffer{rw=true, size=self.volume * self.dim * self.eqn.numWaves * realSize}
	self.rTildeBuf = ctx:buffer{rw=true, size=self.volume * self.dim * self.eqn.numWaves * realSize}
	self.fluxBuf = ctx:buffer{rw=true, size=self.volume * self.dim * self.eqn.numStates * realSize}
	self.derivBuf = ctx:buffer{rw=true, size=self.volume * self.eqn.numStates * realSize}
	
	-- debug only
	self.fluxMatrixBuf = ctx:buffer{rw=true, size=self.volume * self.dim * self.eqn.numStates * self.eqn.numStates * realSize}
	self.orthoErrorBuf = ctx:buffer{rw=true, size=self.volume * self.dim * realSize}
	self.fluxErrorBuf = ctx:buffer{rw=true, size=self.volume * self.dim * realSize}

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

	local ImageGL = require 'cl.imagegl'
	self.texCLMem = ImageGL{context=ctx, tex=self.tex, write=true}
end

function Solver:createCodePrefix()
	local lines = table()
	if dim == 3 then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable'
	end
	if self.app.real == 'double' then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_fp64 : enable'
	end
	
	lines:append(xs:map(function(name,i)
		return '#define gridSize_'..name..' '..tonumber(self.gridSize[name])
	end))
	
	lines:append{
		'#define dim 1',
		'#define numGhost '..self.numGhost,
		'#define xmin '..clnumber(self.mins[1]),
		'#define xmax '..clnumber(self.maxs[1]),
		'#define ymin '..clnumber(self.mins[2]),
		'#define ymax '..clnumber(self.maxs[2]),
		'#define zmin '..clnumber(self.mins[3]),
		'#define zmax '..clnumber(self.maxs[3]),
		'#define dx ((xmax-xmin)/(real)gridSize_x)',
		'#define INDEX(a,b,c)	((a) + gridSize_x * ((b) + gridSize_y * (c)))',
		'#define INDEXV(i)		INDEX((i).x, (i).y, (i).z)',
	}
	lines:append(table{'',2,4,8}:map(function(n)
		return 'typedef '..self.app.real..n..' real'..n..';'
	end))

	lines:append{
		'#define numStates '..self.eqn.numStates,
		'#define numWaves '..self.eqn.numWaves,
		'#define numEigen '..self.eqn.numEigen,
		'constant int4 gridSize = (int4)(gridSize_x, gridSize_y, gridSize_z, 0);',
		'constant real4 dxs = (real4)(dx, 0, 0, 0);',
		'constant int4 stepsize = (int4)(1, gridSize_x, gridSize_x * gridSize_y, gridSize_x * gridSize_y * gridSize_z);',
	}

	if self.eqn.getTypeCode then lines:insert(self.eqn:getTypeCode()) end

	-- run here for teh code, and in buffer for the sizeof()
	lines:insert(self.eqn:getEigenTypeCode())

	-- define i, index, and bounds-check
	lines:insert'#define SETBOUNDS(lhs,rhs)	\\'
	lines:insert'int4 i = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0); \\'
	lines:insert'if (i.x < lhs || i.x >= gridSize_x - rhs \\'
	if self.dim > 1 then lines:insert('|| i.y < lhs || i.y >= gridSize_y - rhs \\') end
	if self.dim > 2 then lines:insert('|| i.z < lhs || i.z >= gridSize_z - rhs \\') end
	lines:insert') return; \\'
	lines:insert'int index = INDEXV(i);'
	

	lines:append(self.displayVars:map(function(var,i)
		return '#define display_'..var.name..' '..i
	end))

	self.codePrefix = lines:concat'\n'
end

function Solver:resetState()
	self.app.cmds:finish()
	self.app.cmds:enqueueNDRangeKernel{kernel=self.initStateKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.app.cmds:finish()
end

-- depends on buffers
function Solver:refreshSolverProgram()

	-- depend on real and nothing else
	-- TODO move to app

	local commonCode = table{
		'typedef '..self.app.real..' real;',
		
		'#define reduce_accum_init INFINITY',
		'#define reduce_operation(x,y) min(x,y)',
		'#define reduce_name reduceMin',
		'#include "reduce.cl"',
		
		'#define reduce_accum_init -INFINITY',
		'#define reduce_operation(x,y) max(x,y)',
		'#define reduce_name reduceMax',
		'#include "reduce.cl"',
	
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

	
	local slopeLimiterCode = 'real slopeLimiter(real r) {'
		.. self.app.slopeLimiters[1+self.slopeLimiter[0]].code 
		.. '}'
	
	local code = table{
		self.codePrefix,
		slopeLimiterCode,
		self.eqn:getEigenCode() or '',
		
		'#define initState_'..self.eqn.initStates[self.initState[0]+1],
		self.eqn:solverCode(clnumber) or '',
		
		'#define calcDisplayVar_dstImage_t '..(self.dim == 3 and 'image3d_t' or 'image2d_t'),
		'#define calcDisplayVar_writeImageArgs '..(dim == 3 and '(int4)(i.x, i.y, i.z, 0)' or '(int2)(i.x, i.y)'),
	
		'#define calcDisplayVar_name calcDisplayVarToTex',
		'#define calcDisplayVar_output_tex',
		'#include "calcDisplayVar.cl"',
		'#undef calcDisplayVar_name',
		'#undef calcDisplayVar_output_tex',
		
		'#define calcDisplayVar_name calcDisplayVarToBuffer',
		'#define calcDisplayVar_output_buffer',
		'#include "calcDisplayVar.cl"',
		'#undef calcDisplayVar_name',
		'#undef calcDisplayVar_output_buffer',

		'#include "solver.cl"',
	}:concat'\n'

	self.solverProgram = require 'cl.program'{context=self.app.ctx, devices={self.app.device}, code=code}

	self.initStateKernel = self.solverProgram:kernel('initState', self.UBuf)

	self.calcDTKernel = self.solverProgram:kernel('calcDT', self.reduceBuf, self.UBuf);
	
	self.calcEigenBasisKernel = self.solverProgram:kernel('calcEigenBasis', self.waveBuf, self.eigenBuf, self.fluxMatrixBuf, self.UBuf)
	self.calcDeltaUTildeKernel = self.solverProgram:kernel('calcDeltaUTilde', self.deltaUTildeBuf, self.UBuf, self.eigenBuf)
	self.calcRTildeKernel = self.solverProgram:kernel('calcRTilde', self.rTildeBuf, self.deltaUTildeBuf, self.waveBuf)
	self.calcFluxKernel = self.solverProgram:kernel('calcFlux', self.fluxBuf, self.UBuf, self.waveBuf, self.eigenBuf, self.deltaUTildeBuf, self.rTildeBuf)
	self.calcDerivFromFluxKernel = self.solverProgram:kernel('calcDerivFromFlux', self.derivBuf, self.fluxBuf)

	self.calcDisplayVarToTexKernel = self.solverProgram:kernel('calcDisplayVarToTex', self.texCLMem)
	self.calcDisplayVarToBufferKernel = self.solverProgram:kernel('calcDisplayVarToBuffer', self.reduceBuf)

	self.calcErrorsKernel = self.solverProgram:kernel('calcErrors', self.orthoErrorBuf, self.fluxErrorBuf, self.waveBuf, self.eigenBuf, self.fluxMatrixBuf)
end

function Solver:refreshBoundaryProgram()
	local lines = table()
	lines:insert(self.codePrefix)
	lines:insert[[
__kernel void boundary(
	__global cons_t* UBuf
) {
	int i = get_global_id(0);
	if (i >= numGhost) return;
]]
	lines:insert(({
		periodic = [[
	UBuf[i] = UBuf[gridSize_x-2*numGhost+i];
]],
		mirror = [[
	UBuf[i] = UBuf[2*numGhost-1-i];
	UBuf[i].mx = -UBuf[i].mx;
]],
		freeflow = [[
	UBuf[i] = UBuf[numGhost];
]],
	})[self.app.boundaryMethods[1+self.boundaryMethods.xmin[0]]])

	lines:insert(({
		periodic = [[
	UBuf[gridSize_x-numGhost+i] = UBuf[numGhost+i];
]],
		mirror = [[
	UBuf[gridSize_x-numGhost+i] = UBuf[gridSize_x-numGhost-1-i];
	UBuf[gridSize_x-numGhost+i].mx = -UBuf[gridSize_x-numGhost+i].mx; 
]],
		freeflow = [[
	UBuf[gridSize_x-numGhost+i] = UBuf[gridSize_x-numGhost-1];
]],
	})[self.app.boundaryMethods[1+self.boundaryMethods.xmax[0]]])
	lines:insert'}'

	local code = lines:concat'\n'
	
	self.app.cmds:finish()
	
	self.boundaryProgram = require 'cl.program'{context=self.app.ctx, devices={self.app.device}, code=code} 
	self.boundaryKernel = self.boundaryProgram:kernel('boundary', self.UBuf);
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
	while reduceSize > 1 do
		iter = iter + 1
		--TODO instead of >> 4, make sure it matches whatever localSize1d is
		-- ... which just so happens to be 16 (i.e. 1 << 4) at the moment
		local nextSize = bit.rshift(reduceSize, 4)
		if 0 ~= bit.band(reduceSize, bit.lshift(1, 4) - 1) then 
			nextSize = nextSize + 1 
		end
		local reduceGlobalSize = math.max(reduceSize, self.localSize1d)
		kernel:setArg(0, src)
		kernel:setArg(2, ffi.new('int[1]', reduceSize))
		kernel:setArg(3, dst)
		self.app.cmds:enqueueNDRangeKernel{kernel=kernel, dim=1, globalSize=reduceGlobalSize, localSize=math.min(reduceGlobalSize, self.localSize1d)}
		self.app.cmds:finish()
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
	self.calcDerivFromFluxKernel:setArg(0, derivBuf)
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDerivFromFluxKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
end

function Solver:integrate(derivBuf, dt)
	-- forward Euler
	self.multAddKernel:setArgs(self.UBuf, self.UBuf, derivBuf, ffi.new('real[1]', dt))
	self.app.cmds:enqueueNDRangeKernel{kernel=self.multAddKernel, globalSize=self.volume * self.eqn.numStates, localSize=self.localSize1d}
end

function Solver:boundary()
	-- 1D:
	self.app.cmds:enqueueNDRangeKernel{kernel=self.boundaryKernel, globalSize=self.localSize1d, localSize=self.localSize1d}
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

	-- integrate flux to state by dt
	self:calcDeriv(self.derivBuf, dt)
	self:integrate(self.derivBuf, dt)

	self.t = self.t + dt
end

function Solver:calcDisplayVarToTex(varIndex)
	local var = self.displayVars[varIndex]
	-- copy to GL
	gl.glFinish()
	self.app.cmds:enqueueAcquireGLObjects{objs={self.texCLMem}}
	self.calcDisplayVarToTexKernel:setArg(1, ffi.new('int[1]', varIndex))
	self.calcDisplayVarToTexKernel:setArg(2, self[var.buffer..'Buf'])
	self.app.cmds:enqueueNDRangeKernel{kernel=self.calcDisplayVarToTexKernel, dim=self.dim, globalSize=self.gridSize:ptr(), localSize=self.localSize:ptr()}
	self.app.cmds:enqueueReleaseGLObjects{objs={self.texCLMem}}
	self.app.cmds:finish()
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

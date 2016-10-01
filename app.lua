#!/usr/bin/env luajit
local class = require 'ext.class'
local ImGuiApp = require 'imguiapp'
local ig = require 'ffi.imgui'
local gl = require 'ffi.OpenGL'
local cl = require 'ffi.OpenCL'
local ffi = require 'ffi'
local bit = require 'bit'
local glreport = require 'gl.report'

local vec3sz = require 'ffi.vec.create_ffi'(3,'size_t','sz')

local gridSize = vec3sz(256,1,1)
local volume = tonumber(gridSize:volume())
local dim = 1
local cfl = .5
local numStates = 3
local numWaves = 3
local numEigen = 2 * 3 * 3

local xmin = -1
local xmax = 1
local gamma = 7/5

local offset = vec3sz(0,0,0)
local localSize1d = 16
local localSize = dim < 3 and vec3sz(16,16,16) or vec3sz(8,8,8) 

local real = 'float'
local realSize = ffi.sizeof(real)
ffi.cdef('typedef '..real..' real;')

local primType = 'prim_t'
local primTypeCode = [[
typedef struct {
	real rho;
	real vx;
	real P;
} ]]..primType..';\n'
ffi.cdef(primTypeCode)
local primSize = ffi.sizeof(primType)

local consType = 'cons_t'
local consTypeCode = [[
typedef struct {
	real rho;
	real mx;
	real ETotal;
} ]]..consType..';\n'
ffi.cdef(consTypeCode)
local consSize = ffi.sizeof(consType)
assert(realSize * numStates == consSize)

local displayVars = table{'rho', 'vx', 'P', 'eInt', 'eKin', 'eTotal'}

local function clnumber(x)
	local s = tostring(x)
	if s:find'e' then return s end
	if not s:find('%.') then s = s .. '.' end
	return s
end

local HydroCLApp = class(ImGuiApp)

function HydroCLApp:initGL(...)
	HydroCLApp.super.initGL(self, ...)

	-- TODO generate the codeprogram
	local lines = table()
	if dim == 3 then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable'
	end
	if real == 'double' then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_fp64 : enable'
	end
	
	lines:append(table{'x','y','z'}:map(function(name,i)
		return '#define gridSize_'..name..' '..gridSize[name]
	end))
	
	lines:append{
		'#define dim 1',
		'#define numStates '..numStates,
		'#define numWaves '..numWaves,
		'#define numEigen '..numEigen,
		'#define xmin '..clnumber(xmin),
		'#define xmax '..clnumber(xmax),
		'#define dx ((xmax-xmin)/(real)gridSize_x)',
		'#define INDEX(a,b,c)	((a) + gridSize_x * ((b) + gridSize_y * (c)))',
		'#define INDEXV(i)		INDEX((i).x, (i).y, (i).z)',
		'#define gamma '..clnumber(gamma),
		'#define WRITEIMAGEARGS '..(dim == 3 and '(int4)(i.x, i.y, i.z, 0)' or '(int2)(i.x, i.y)'),
		'typedef '..(dim == 3 and 'image3d_t' or 'image2d_t')..' dstimage_t;',
	}
	lines:append(table{'',2,4,8}:map(function(n)
		return 'typedef '..real..n..' real'..n..';'
	end))
	lines:insert(primTypeCode)
	lines:insert(consTypeCode)

	-- define i, index, and bounds-check
	lines:insert'#define SETBOUNDS(lhs,rhs)	\\'
	lines:insert'int4 i = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0); \\'
	lines:insert'if (i.x < lhs || i.x >= gridSize_x - rhs \\'
	if dim > 1 then lines:insert('|| i.y < lhs || i.y >= gridSize_y - rhs \\') end
	if dim > 2 then lines:insert('|| i.z < lhs || i.z >= gridSize_z - rhs \\') end
	lines:insert') return; \\'
	lines:insert'int index = INDEXV(i);'
	lines:append(displayVars:map(function(var,i)
		return '#define display_'..var..' '..i
	end))
	lines:insert'#include "solver.cl"'	
	
	local code = lines:concat'\n'

	self.platform, self.device, self.ctx, self.cmds, self.program
	= require 'cl'{
		device={gpu=true}, 
		context={glSharing=true}, 
		program={code=code},
	}

	self.UBuf = self.ctx:buffer{rw=true, size=volume*consSize}
	self.dtBuf = self.ctx:buffer{rw=true, size=volume*realSize}
	self.reduceFinalMem = ffi.new(real..'[1]', 0)
	self.dtSwapBuf = self.ctx:buffer{rw=true, size=volume*realSize/localSize1d}
	self.waveBuf = self.ctx:buffer{rw=true, size=volume*dim*numWaves*realSize}
	self.eigenBuf = self.ctx:buffer{rw=true, size=volume*dim*numEigen*realSize}
	self.deltaUTildeBuf = self.ctx:buffer{rw=true, size=volume*dim*numWaves*realSize}
	self.rTildeBuf = self.ctx:buffer{rw=true, size=volume*dim*numWaves*realSize}
	self.fluxBuf = self.ctx:buffer{rw=true, size=volume*dim*consSize}
	self.derivBuf = self.ctx:buffer{rw=true, size=volume*consSize}

	self.calcDTKernel = self.program:kernel('calcDT', self.dtBuf, self.UBuf);
	self.reduceMinKernel = self.program:kernel('reduceMin',
		self.dtBuf, 
		-- hmm, how to handle local space declarations?
		{ptr=nil, size=localSize1d*realSize},
		ffi.new('int[1]', volume),
		self.dtSwapBuf)

	self.calcEigenBasisKernel = self.program:kernel('calcEigenBasis', self.waveBuf, self.eigenBuf, self.UBuf)
	self.calcDeltaUTildeKernel = self.program:kernel('calcDeltaUTilde', self.deltaUTildeBuf, self.UBuf, self.eigenBuf)
	self.calcRTildeKernel = self.program:kernel('calcRTilde', self.rTildeBuf, self.deltaUTildeBuf, self.waveBuf)
	self.calcFluxKernel = self.program:kernel('calcFlux', self.fluxBuf, self.UBuf, self.waveBuf, self.eigenBuf, self.deltaUTildeBuf, self.rTildeBuf)
	self.calcDerivFromFluxKernel = self.program:kernel('calcDerivFromFlux', self.derivBuf, self.fluxBuf)
	self.multAddToKernel = self.program:kernel'multAddTo'

	local initStateKernel = self.program:kernel('initState', self.UBuf)
	self.cmds:enqueueNDRangeKernel{kernel=initStateKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.cmds:finish()
	
	local GLTex2D = require 'gl.tex2d'
	self.tex = GLTex2D{
		width = gridSize.x,
		height = 1,
		internalFormat = gl.GL_RGBA32F,
		format = gl.GL_RGBA,
		type = gl.GL_FLOAT,
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_LINEAR,
		wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT},
	}

	local ImageGL = require 'cl.imagegl'
	self.texCLMem = ImageGL{context=self.ctx, tex=self.tex, write=true}

	self.convertToTexKernel = self.program:kernel('convertToTex', self.texCLMem, self.UBuf)

	local GLProgram = require 'gl.program'
	local graphShaderCode = file['graph.shader']
	self.graphShader = GLProgram{
		vertexCode = '#define VERTEX_SHADER\n'..graphShaderCode,
		fragmentCode = '#define FRAGMENT_SHADER\n'..graphShaderCode,
		uniforms = {
			'tex', 
			'scale', 
			'xmin', 'xmax', 
			'useLog', 
			'axis', 
			'ambient',
			'size', 
		},
	}
	self.graphShader:use()
	gl.glUniform1i(self.graphShader.uniforms.tex, 0)
	gl.glUniform1f(self.graphShader.uniforms.scale, 1)
	gl.glUniform2f(self.graphShader.uniforms.xmin, xmin, 0)
	gl.glUniform2f(self.graphShader.uniforms.xmax, xmax, 0)
	gl.glUniform1i(self.graphShader.uniforms.useLog, false)
	gl.glUniform1i(self.graphShader.uniforms.axis, dim)
	gl.glUniform1f(self.graphShader.uniforms.ambient, 1)
	gl.glUniform2f(self.graphShader.uniforms.size, gridSize.x, gridSize.y)
	self.graphShader:useNone()
end

function HydroCLApp:reduceMinDT()
	local reduceSize = volume
	local dst = self.dtSwapBuf
	local src = self.dtBuf
	while reduceSize > 1 do
		--TODO instead of >> 4, make sure it matches whatever localSize1d is
		-- ... which just so happens to be 16 (i.e. 1 << 4) at the moment
		local nextSize = bit.rshift(reduceSize, 4)
		if 0 ~= bit.band(reduceSize, bit.lshift(1, 4) - 1) then 
			nextSize = nextSize + 1 
		end
		local reduceGlobalSize = math.max(reduceSize, localSize1d)
		self.reduceMinKernel:setArg(0, src)
		self.reduceMinKernel:setArg(2, ffi.new('int[1]',reduceSize))
		self.reduceMinKernel:setArg(3, dst)
		self.cmds:enqueueNDRangeKernel{kernel=self.reduceMinKernel, dim=1, globalSize=reduceGlobalSize, localSize=math.min(reduceGlobalSize, localSize1d)}
		self.cmds:finish()
		dst, src = src, dst
		reduceSize = nextSize
	end
	self.cmds:enqueueReadBuffer{buffer=src, block=true, size=realSize, ptr=self.reduceFinalMem}
	return self.reduceFinalMem[0]
end

function HydroCLApp:calcDeriv(derivBuf, dt)
	self.cmds:enqueueNDRangeKernel{kernel=self.calcEigenBasisKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.cmds:enqueueNDRangeKernel{kernel=self.calcDeltaUTildeKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.cmds:enqueueNDRangeKernel{kernel=self.calcRTildeKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.calcFluxKernel:setArg(6, ffi.new('real[1]', dt))
	self.cmds:enqueueNDRangeKernel{kernel=self.calcFluxKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.calcDerivFromFluxKernel:setArg(0, derivBuf)
	self.cmds:enqueueNDRangeKernel{kernel=self.calcDerivFromFluxKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
end

function HydroCLApp:integrate(x, dx_dt, dt)
	-- forward Euler
	self.multAddToKernel:setArgs(x, dx_dt, ffi.new('real[1]', dt))
	self.cmds:enqueueNDRangeKernel{kernel=self.multAddToKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
end

local vec3 = require 'vec.vec3'
local displayVarInfo = displayVars:map(function(var,i)
	return {
		enabled = ffi.new('bool[1]', i==1),
		color = vec3(math.random(), math.random(), math.random()):normalize(),
	}
end)

function HydroCLApp:update(...)

	-- calc cell wavespeeds -> dts
	local dt = self.fixedDT
	if not dt then
		self.cmds:enqueueNDRangeKernel{kernel=self.calcDTKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
		dt = cfl * self:reduceMinDT()
	end

	-- integrate flux to state by dt
	self:calcDeriv(self.derivBuf, dt)
	self:integrate(self.UBuf, self.derivBuf, dt)

	gl.glClearColor(0,0,0,1)
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	
	local ar = self.width / self.height
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(-ar, ar, -1, 1, -1, 1)

	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
	gl.glScalef(1, 1/2, 1)

	for i,info in ipairs(displayVarInfo) do
		if info.enabled[0] then
			self:graphVar(i, info.color)
		end
	end

	HydroCLApp.super.update(self, ...)
end

function HydroCLApp:graphVar(i, color)
	-- copy to GL
	gl.glFinish()
	self.cmds:enqueueAcquireGLObjects{objs={self.texCLMem}}
	self.convertToTexKernel:setArg(2, ffi.new('int[1]', i))
	self.cmds:enqueueNDRangeKernel{kernel=self.convertToTexKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.cmds:enqueueReleaseGLObjects{objs={self.texCLMem}}
	self.cmds:finish()

	-- display

	self.graphShader:use()
	self.tex:bind()

	gl.glColor3f(table.unpack(color))
	gl.glBegin(gl.GL_LINE_STRIP)
	local step = 1
	for i=2,tonumber(gridSize.x)-2,step do
		local x = (i+.5)/tonumber(gridSize.x)
		gl.glVertex2f(x, 0)
	end
	gl.glEnd()
	--[[
	gl.glBegin(gl.GL_QUADS)
	for k,v in ipairs{{0,0}, {1,0}, {1,1}, {0,1}} do
		gl.glTexCoord2f(table.unpack(v))
		gl.glVertex2f(v[1]*2-1, v[2]*2-1)
	end
	gl.glEnd()
	--]]
	
	self.tex:unbind()
	self.graphShader:useNone()
end

function HydroCLApp:updateGUI()
	ig.igText('variables:')
	for i,var in ipairs(displayVars) do
		ig.igCheckbox(var, displayVarInfo[i].enabled)
	end
end

return HydroCLApp

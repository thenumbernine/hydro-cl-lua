#!/usr/bin/env luajit
local class = require 'ext.class'
local ImGuiApp = require 'imguiapp'
local ffi = require 'ffi'
local ig = require 'ffi.imgui'
local gl = require 'ffi.OpenGL'
local cl = require 'ffi.OpenCL'
local sdl = require 'ffi.sdl'
local bit = require 'bit'
local glreport = require 'gl.report'

local vec3sz = require 'ffi.vec.create_ffi'(3,'size_t','sz')

local gridSize = vec3sz(256,1,1)
local volume = tonumber(gridSize:volume())
local dim = 1
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

local function clnumber(x)
	local s = tostring(x)
	if s:find'e' then return s end
	if not s:find('%.') then s = s .. '.' end
	return s
end

local HydroCLApp = class(ImGuiApp)

HydroCLApp.title = 'Hydrodynamics in OpenCL'

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
	

	self.displayVars = table()
	local vec3 = require 'vec.vec3'
	local function makevars(buffer, ...)
		for i=1,select('#',...) do
			local name = tostring(select(i,...))
			self.displayVars:insert{
				buffer = buffer,
				name = buffer..'_'..name,
				enabled = ffi.new('bool[1]', buffer == 'U' and i <= 2),
				color = vec3(math.random(), math.random(), math.random()):normalize(),
			}
		end
	end
	makevars('U', 'rho', 'vx', 'P', 'eInt', 'eKin', 'eTotal')
	makevars('wave', range(0,numWaves-1):unpack())
	makevars('eigen', range(0,numEigen-1):unpack())
	makevars('dt', '0')
	makevars('deltaUTilde', range(0,numWaves-1):unpack())
	makevars('rTilde', range(0,numWaves-1):unpack())
	makevars('flux', range(0,numStates-1):unpack())
	makevars('deriv', range(0,numStates-1):unpack())
	makevars('orthoError', '0')
	makevars('fluxError', '0')

	lines:append(self.displayVars:map(function(var,i)
		return '#define display_'..var.name..' '..i
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
	self.reduceFinalMem = ffi.new('real[1]', 0)
	self.dtSwapBuf = self.ctx:buffer{rw=true, size=volume*realSize/localSize1d}
	self.waveBuf = self.ctx:buffer{rw=true, size=volume*dim*numWaves*realSize}
	self.eigenBuf = self.ctx:buffer{rw=true, size=volume*dim*numEigen*realSize}
	self.deltaUTildeBuf = self.ctx:buffer{rw=true, size=volume*dim*numWaves*realSize}
	self.rTildeBuf = self.ctx:buffer{rw=true, size=volume*dim*numWaves*realSize}
	self.fluxBuf = self.ctx:buffer{rw=true, size=volume*dim*consSize}
	self.derivBuf = self.ctx:buffer{rw=true, size=volume*consSize}
	
	-- debug only
	self.fluxMatrixBuf = self.ctx:buffer{rw=true, size=volume*dim*numStates*numStates*realSize}
	self.orthoErrorBuf = self.ctx:buffer{rw=true, size=volume*dim*realSize}
	self.fluxErrorBuf = self.ctx:buffer{rw=true, size=volume*dim*realSize}

	self.calcDTKernel = self.program:kernel('calcDT', self.dtBuf, self.UBuf);
	self.reduceMinKernel = self.program:kernel('reduceMin',
		self.dtBuf, 
		-- hmm, how to handle local space declarations?
		{ptr=nil, size=localSize1d*realSize},
		ffi.new('int[1]', volume),
		self.dtSwapBuf)

	self.calcEigenBasisKernel = self.program:kernel('calcEigenBasis', self.waveBuf, self.eigenBuf, self.fluxMatrixBuf, self.UBuf)
	self.calcDeltaUTildeKernel = self.program:kernel('calcDeltaUTilde', self.deltaUTildeBuf, self.UBuf, self.eigenBuf)
	self.calcRTildeKernel = self.program:kernel('calcRTilde', self.rTildeBuf, self.deltaUTildeBuf, self.waveBuf)
	self.calcFluxKernel = self.program:kernel('calcFlux', self.fluxBuf, self.UBuf, self.waveBuf, self.eigenBuf, self.deltaUTildeBuf, self.rTildeBuf)
	self.calcDerivFromFluxKernel = self.program:kernel('calcDerivFromFlux', self.derivBuf, self.fluxBuf)
	self.multAddKernel = self.program:kernel'multAdd'

	self.calcErrorsKernel = self.program:kernel('calcErrors', self.orthoErrorBuf, self.fluxErrorBuf, self.waveBuf, self.eigenBuf, self.fluxMatrixBuf)

	self.initStateKernel = self.program:kernel('initState', self.UBuf)
	self:resetState()
	
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

	self.convertToTexKernel = self.program:kernel('convertToTex', self.texCLMem)

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

function HydroCLApp:resetState()
	self.cmds:finish()
	self.cmds:enqueueNDRangeKernel{kernel=self.initStateKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.cmds:finish()
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
	self.cmds:enqueueNDRangeKernel{kernel=self.calcErrorsKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.cmds:enqueueNDRangeKernel{kernel=self.calcDeltaUTildeKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.cmds:enqueueNDRangeKernel{kernel=self.calcRTildeKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.calcFluxKernel:setArg(6, ffi.new('real[1]', dt))
	self.cmds:enqueueNDRangeKernel{kernel=self.calcFluxKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.calcDerivFromFluxKernel:setArg(0, derivBuf)
	self.cmds:enqueueNDRangeKernel{kernel=self.calcDerivFromFluxKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
end

function HydroCLApp:integrate(derivBuf, dt)
	-- forward Euler
	self.multAddKernel:setArgs(self.UBuf, self.UBuf, derivBuf, ffi.new('real[1]', dt))
	self.cmds:enqueueNDRangeKernel{kernel=self.multAddKernel, globalSize=volume*numStates, localSize=localSize1d}
end

local xScale = ffi.new('float[1]', 1)
local yScale = ffi.new('float[1]', .5)

local updateMethod
local useFixedDT = ffi.new('bool[1]', false)
local fixedDT = ffi.new('float[1]', 0)
local currentDT = ffi.new('float[1]', 0)
local cfl = ffi.new('float[1]', .5)

function HydroCLApp:update(...)

	if updateMethod == 'reset' then
		self:resetState()
		updateMethod = nil 
	end

	if updateMethod then
		if updateMethod == 'step' then 
			print('performing single step...')
			updateMethod = nil 
		end

		-- calc cell wavespeeds -> dts
		if useFixedDT[0] then
			dt = tonumber(fixedDT[0])
		else
			self.cmds:enqueueNDRangeKernel{kernel=self.calcDTKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
			dt = tonumber(cfl[0]) * self:reduceMinDT()
		end
		currentDT[0] = dt
		
		-- integrate flux to state by dt
		self:calcDeriv(self.derivBuf, dt)
		self:integrate(self.derivBuf, dt)
	end

	gl.glClearColor(.3,.2,.5,1)
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	
	local ar = self.width / self.height
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(-ar, ar, -1, 1, -1, 1)

	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
	gl.glScalef(xScale[0], yScale[0], 1)

	for i,var in ipairs(self.displayVars) do
		if var.enabled[0] then
			self:renderDisplayVar(i, var)
		end
	end

	HydroCLApp.super.update(self, ...)
end

function HydroCLApp:renderDisplayVar(i, var)
	-- copy to GL
	gl.glFinish()
	self.cmds:enqueueAcquireGLObjects{objs={self.texCLMem}}
	self.convertToTexKernel:setArg(1, ffi.new('int[1]', i))
	self.convertToTexKernel:setArg(2, self[var.buffer..'Buf'])
	self.cmds:enqueueNDRangeKernel{kernel=self.convertToTexKernel, dim=dim, globalSize=gridSize:ptr(), localSize=localSize:ptr()}
	self.cmds:enqueueReleaseGLObjects{objs={self.texCLMem}}
	self.cmds:finish()

	-- display

	self.graphShader:use()
	self.tex:bind()

	gl.glColor3f(table.unpack(var.color))
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
	if ig.igButton(updateMethod and 'Stop' or 'Start') then
		updateMethod = not updateMethod
	end
	if ig.igButton'Step' then
		updateMethod = 'step'
	end
	
	ig.igCheckbox('use fixed dt', useFixedDT)
	ig.igInputFloat('fixed dt', fixedDT)
	ig.igInputFloat('current dt', currentDT)
	ig.igInputFloat('cfl', cfl)

	ig.igSliderFloat('x scale', xScale, 0, 100, '%.3f', 10)
	ig.igSliderFloat('y scale', yScale, 0, 100, '%.3f', 10)

	if ig.igCollapsingHeader'variables:' then
		local lastSection
		local sectionEnabled
		for i,var in ipairs(self.displayVars) do
			local section = var.buffer
			if section ~= lastSection then
				sectionEnabled = ig.igCollapsingHeader(section..' variables:')
			end
			if sectionEnabled then
				ig.igCheckbox(var.name, var.enabled)
			end
			lastSection = section
		end
	end
end

function HydroCLApp:event(event, ...)
	HydroCLApp.super.event(self, event, ...)
	if ig.igGetIO()[0].WantCaptureKeyboard then return end
	if event.type == sdl.SDL_KEYDOWN then
		if event.key.keysym.sym == sdl.SDLK_SPACE then
			updateMethod = not updateMethod
		elseif event.key.keysym.sym == ('u'):byte() then
			updateMethod = 'step'
		elseif event.key.keysym.sym == ('r'):byte() then
			print'resetting...'
			updateMethod = 'reset'
		end
	end
end

return HydroCLApp

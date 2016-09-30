#!/usr/bin/env luajit
local class = require 'ext.class'
local ImGuiApp = require 'imguiapp'
local ig = require 'ffi.imgui'
local gl = require 'ffi.OpenGL'
local cl = require 'ffi.OpenCL'
local ffi = require 'ffi'

local glreport = require 'gl.report'

local HydroCLApp = class(ImGuiApp)

local gridsize = 256
local dim = 1
local numStates = 3
local real = 'float'
local realsize = ffi.sizeof(real)

local xmin = -1
local xmax = 1
local gamma = 7/5

local function clnumber(x)
	local s = tostring(x)
	if s:find'e' then return s end
	if not s:find('%.') then s = s .. '.' end
	return s
end

function HydroCLApp:initGL(...)
	glreport'here'
	HydroCLApp.super.initGL(self, ...)
	glreport'here'


	-- TODO generate the codeprogram
	local lines = table()
	if dim == 3 then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable'
	end
	if real == 'double' then
		lines:insert'#pragma OPENCL EXTENSION cl_khr_fp64 : enable'
	end
	
	lines:append{
		'#define SIZE_X '..gridsize,
		'#define SIZE_Y 1',
		'#define SIZE_Z 1',
		'#define DIM 1',
		'#define real '..real,
		'#define numStates '..numStates,
		'#define xmin '..clnumber(xmin),
		'#define xmax '..clnumber(xmax),
		'#define INDEX(a,b,c)	((a) + SIZE_X * ((b) + SIZE_Y * (c)))',
		'#define INDEXV(i)		INDEX((i).x, (i).y, (i).z)',
		'#define gamma '..clnumber(gamma),
		'#define IMAGETYPE '..(dim == 3 and 'image3d_t' or 'image2d_t'),
		'#define WRITEIMAGEARGS '..(dim == 3 and '(int4)(i.x, i.y, i.z, 0)' or '(int2)(i.x, i.y)'),
	}

	-- define i, index, and bounds-check
	lines:insert'#define BEGIN_KERNEL()	\\'
	lines:insert'int4 i = (int4)(get_global_id(0), get_global_id(1), get_global_id(2), 0); \\'
	lines:insert'if (i.x < 0 || i.x >= SIZE_X \\'
	if dim > 1 then lines:insert('|| i.y < 0 || i.y >= SIZE_Y \\') end
	if dim > 2 then lines:insert('|| i.z < 0 || i.z >= SIZE_Z \\') end
	lines:insert') return; \\'
	lines:insert'int index = INDEXV(i);'
	
	local code = lines:concat'\n'..'\n'..file['solver.cl']

	self.platform, self.device, self.ctx, self.cmds, self.program
	= require 'cl'{
		device={gpu=true}, 
		context={glSharing=true}, 
		program={code=code},
	}

	-- conservative/state buffer
	self.usBuf = self.ctx:buffer{rw=true, size=gridsize*numStates*realsize}
	self.fsBuf = self.ctx:buffer{rw=true, size=gridsize*numStates*realsize}

	local initStateKernel = self.program:kernel('initState', self.usBuf)
	self.cmds:enqueueNDRangeKernel{kernel=initStateKernel, globalSize={gridsize}, localSize={16}}
	self.cmds:finish()

	glreport'here'
	self.tex = require 'gl.tex2d'{
		width = gridsize,
		height = 1,
		internalFormat = gl.GL_RGBA32F,
		format = gl.GL_RGBA,
		type = gl.GL_FLOAT,
		minFilter = gl.GL_NEAREST,
		magFilter = gl.GL_LINEAR,
		wrap = {s=gl.GL_REPEAT, t=gl.GL_REPEAT},
--		data = ffi.new('float[?]', gridsize * 1 * 4 * ffi.sizeof'float'),
	}
	glreport'here'

	self.hsvTex= require 'gl.hsvtex'(1024)

	local ImageGL = require 'cl.imagegl'
	self.texCLMem = ImageGL{context=self.ctx, tex=self.tex, write=true}

	self.convertToTexKernel = self.program:kernel('convertToTex', self.texCLMem, self.usBuf)
	
end

function HydroCLApp:update(...)

	-- calc cell wavespeeds -> dts
	-- reduce to find min dt, scale by cfl
	-- calc flux
	-- integrate flux to state by dt

--	self.cmds:enqueueNDRangeKernel{kernel=self.calcFluxKernel, globalSize={gridsize}, localSize={16}}
--	self.cmds:finish()

	-- copy to GL
	gl.glFinish()
	self.cmds:enqueueAcquireGLObjects{objs={self.texCLMem}}
	self.cmds:enqueueNDRangeKernel{kernel=self.convertToTexKernel, globalSize={gridsize}, localSize={16}}
	self.cmds:enqueueReleaseGLObjects{objs={self.texCLMem}}
	self.cmds:finish()
	glreport'here'

	-- display

	gl.glClearColor(.4, .8, .8, 1)
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	
	local ar = self.width / self.height
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(-ar, ar, -1, 1, -1, 1)

	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()

	self.tex:enable()
	self.tex:bind()	
	gl.glBegin(gl.GL_QUADS)
	for k,v in ipairs{{0,0}, {1,0}, {1,1}, {0,1}} do
		gl.glTexCoord2f(table.unpack(v))
		gl.glVertex2f(v[1]*2-1, v[2]*2-1)
	end
	gl.glEnd()
	self.tex:unbind()
	self.tex:disable()

	HydroCLApp.super.update(self, ...)

end

function HydroCLApp:updateGUI()
	ig.igText('Hello, world!')
end

return HydroCLApp

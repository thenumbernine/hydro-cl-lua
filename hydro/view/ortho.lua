local class = require 'ext.class'
local gl = require 'gl'
local vec2d = require 'vec-ffi.vec2d'
local matrix_ffi = require 'matrix.ffi'

matrix_ffi.real = 'float'	-- default matrix_ffi type

local OrthoView = class()

function OrthoView:init()
	local zoom = cmdline.zoom or .5
	self.zoom = vec2d(zoom, zoom)
	self.pos = vec2d()

	self.modelViewMatrix = matrix_ffi.zeros(4,4)
	self.projectionMatrix = matrix_ffi.zeros(4,4)
	self.modelViewProjectionMatrix = matrix_ffi.zeros(4,4)
end

-- returns xmin, xmax, ymin, ymax, zmin, zmax
function OrthoView:getOrthoBounds(ar)
	return 
		self.pos.x - ar * .5 / self.zoom.x, 
		self.pos.x + ar * .5 / self.zoom.x,
		self.pos.y - .5 / self.zoom.y,
		self.pos.y + .5 / self.zoom.y,
		-1, 1
end

function OrthoView:projection(ar)
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(self:getOrthoBounds(ar))
end

function OrthoView:modelview()
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
end

function OrthoView:setup(ar)
	self:projection(ar)
	self:modelview()

	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.modelViewMatrix.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.projectionMatrix.ptr)
	self.modelViewProjectionMatrix:mul(self.projectionMatrix, self.modelViewMatrix)
end

function OrthoView:mousePan(dx, dy, screenWidth, screenHeight)
	local ar = screenWidth / screenHeight
	self.pos = self.pos + vec2d(
		-dx * ar / screenWidth / self.zoom.x,
		dy / screenHeight / self.zoom.y)
end

function OrthoView:mouseZoom(dx, dy)
	self.zoom.x = self.zoom.x * math.exp(-dx * -.03)
	self.zoom.y = self.zoom.y * math.exp(dy * -.03)
end

return OrthoView 

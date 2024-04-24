local class = require 'ext.class'
local gl = require 'gl'
local vec2d = require 'vec-ffi.vec2d'
local matrix_ffi = require 'matrix.ffi'


local OrthoView = class()

function OrthoView:init()
	local zoom = cmdline.zoom or .5
	self.zoom = vec2d(zoom, zoom)
	self.pos = vec2d()

	self.mvMat = matrix_ffi.zeros({4,4}, 'float')
	self.projMat = matrix_ffi.zeros({4,4}, 'float')
	self.mvProjMat = matrix_ffi.zeros({4,4}, 'float')
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

function OrthoView:setupProjection(ar)
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(self:getOrthoBounds(ar))
end

function OrthoView:setupModelView()
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
end

function OrthoView:setup(ar)
	self:setupProjection(ar)
	self:setupModelView()

	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.mvMat.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.projMat.ptr)
	self.mvProjMat:mul(self.projMat, self.mvMat)
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

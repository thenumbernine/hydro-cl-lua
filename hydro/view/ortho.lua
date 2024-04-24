local class = require 'ext.class'
local gl = require 'gl'
local vec2d = require 'vec-ffi.vec2d'
local vec3d = require 'vec-ffi.vec3d'
local quatd = require 'vec-ffi.quatd'

local OrthoView = class()

OrthoView.znear = -1
OrthoView.zfar = 1

function OrthoView:init(args)
	-- in glapp.view, zoom is replaced with orthoSize, a scalar
	self.zoom = vec2d(.5, .5)
	self.pos = vec2d()
	self.orbit = vec3d(0,0,0)	-- orbit center ... not used yet
	self.angle = quatd(0,0,0,1)	-- ortho angle ... not used yet
	if args then
		local zoom = tonumber(args.zoom)
		if zoom then
			self.zoom:set(zoom, zoom)
		end
	end

	local matrix_ffi = require 'matrix.ffi'
	self.mvMat = matrix_ffi.zeros({4,4}, 'float')
	self.projMat = matrix_ffi.zeros({4,4}, 'float')
	self.mvProjMat = matrix_ffi.zeros({4,4}, 'float')
end

-- returns xmin, xmax, ymin, ymax, zmin, zmax
function OrthoView:getOrthoBounds(aspectRatio)
	return 
		self.pos.x - aspectRatio * .5 / self.zoom.x, 
		self.pos.x + aspectRatio * .5 / self.zoom.x,
		self.pos.y - .5 / self.zoom.y,
		self.pos.y + .5 / self.zoom.y,
		self.znear,
		self.zfar
end

function OrthoView:setupProjection(aspectRatio)
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glOrtho(self:getOrthoBounds(aspectRatio))
end

function OrthoView:setupModelView()
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
end

function OrthoView:setup(aspectRatio)
	self:setupProjection(aspectRatio)
	self:setupModelView()

	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.mvMat.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.projMat.ptr)
	self.mvProjMat:mul(self.projMat, self.mvMat)
end

-- not in glapp.view
function OrthoView:mousePan(dx, dy, screenWidth, screenHeight)
	local aspectRatio = screenWidth / screenHeight
	self.pos = self.pos + vec2d(
		-dx * aspectRatio / screenWidth / self.zoom.x,
		dy / screenHeight / self.zoom.y)
end

-- not in glapp.view
function OrthoView:mouseZoom(dx, dy)
	self.zoom.x = self.zoom.x * math.exp(-dx * -.03)
	self.zoom.y = self.zoom.y * math.exp(dy * -.03)
end

return OrthoView 

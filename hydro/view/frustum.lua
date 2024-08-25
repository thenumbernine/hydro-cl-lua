local class = require 'ext.class'
local table = require 'ext.table'
local gl = require 'gl'
local vec3d = require 'vec-ffi.vec3d'
local quatd = require 'vec-ffi.quatd'

local FrustumView = class()

FrustumView.zfar = 1000
FrustumView.znear = .1

function FrustumView:init(args)
	self.pos = vec3d()
	self.dist = 3	-- TODO integrate this with .orbit ...
	self.orbit = vec3d(0,0,0)	-- orbit center ... not used yet
	self.angle = quatd(0,0,0,1)
	if args then
		if args.dist then self.dist = args.dist end
		if args.pos then self.pos:set(unpack(args.pos)) end
		if args.orbit then self.orbit:set(unpack(args.orbit)) end
		if args.angle then self.angle:set(unpack(args.angle)):normalize(self.angle) end
	end

	local matrix_ffi = require 'matrix.ffi'
	self.mvMat = matrix_ffi.zeros({4,4}, 'float')
	self.projMat = matrix_ffi.zeros({4,4}, 'float')
	self.mvProjMat = matrix_ffi.zeros({4,4}, 'float')
end

function FrustumView:getFrustumBounds(aspectRatio)
	return
		-aspectRatio * self.znear,
		aspectRatio * self.znear,
		-self.znear,
		self.znear,
		self.znear,
		self.zfar
end

function FrustumView:setupProjection(aspectRatio)
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glFrustum(self:getFrustumBounds(aspectRatio))
end

function FrustumView:setupModelView()
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
	gl.glTranslatef(0,0,-self.dist)
	local angleAxis = self.angle:toAngleAxis()
	gl.glRotatef(-angleAxis.w, angleAxis.x, angleAxis.y, angleAxis.z)
	gl.glTranslatef(-self.pos.x, -self.pos.y, -self.pos.z)
end

function FrustumView:setup(aspectRatio)
	self:setupProjection(aspectRatio)
	self:setupModelView()

	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.mvMat.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.projMat.ptr)
	-- TODO :mul() is transposed from :mul4x4()
	self.mvProjMat:mul4x4(self.projMat, self.mvMat)
end

-- not in glapp.view
function FrustumView:mousePan(dx, dy, screenWidth, screenHeight)
-- [[ frustum pan
	local magn = math.sqrt(dx * dx + dy * dy)
	local fdx = dx / magn
	local fdy = dy / magn
	local rotation = quatd():fromAngleAxis(-fdy, -fdx, 0, magn)
	self.angle = (self.angle * rotation):normalize()
--]]
--[[ ortho pan
	local aspectRatio = screenWidth / screenHeight
	self.pos = self.pos + vec3d(
		-dx * aspectRatio / screenWidth * self.dist,
		dy / screenHeight * self.dist,
		0)
--]]
end

-- not in glapp.view
function FrustumView:mouseZoom(dx, dy)
	self.dist = self.dist * math.exp(dy * -.03)
end

return FrustumView

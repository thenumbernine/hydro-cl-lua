local class = require 'ext.class'
local table = require 'ext.table'
local gl = require 'gl'
local vec3d = require 'vec-ffi.vec3d'
local quatd = require 'vec-ffi.quatd'
local matrix_ffi = require 'matrix.ffi'


local FrustumView = class()

function FrustumView:init()
	self.dist = cmdline.frustumDist or 3
	self.pos = vec3d()
	self.angle = quatd(table.unpack(cmdline.frustumAngle or {0,0,0,1}))
	self.angle:normalize(self.angle)

	self.modelViewMatrix = matrix_ffi.zeros({4,4}, 'float')
	self.projectionMatrix = matrix_ffi.zeros({4,4}, 'float')
	self.modelViewProjectionMatrix = matrix_ffi.zeros({4,4}, 'float')
end

FrustumView.zFar = 1000
FrustumView.zNear = .1
function FrustumView:projection(ar)
	gl.glMatrixMode(gl.GL_PROJECTION)
	gl.glLoadIdentity()
	gl.glFrustum(
		-ar * self.zNear, 
		ar * self.zNear,
		-self.zNear,
		self.zNear,
		self.zNear,
		self.zFar)
end

function FrustumView:modelview()
	gl.glMatrixMode(gl.GL_MODELVIEW)
	gl.glLoadIdentity()
	gl.glTranslatef(0,0,-self.dist)
	local angleAxis = self.angle:toAngleAxis()
	gl.glRotatef(-angleAxis.w, angleAxis.x, angleAxis.y, angleAxis.z)
	gl.glTranslatef(-self.pos.x, -self.pos.y, -self.pos.z)
end

function FrustumView:setup(ar)
	self:projection(ar)
	self:modelview()

	gl.glGetFloatv(gl.GL_MODELVIEW_MATRIX, self.modelViewMatrix.ptr)
	gl.glGetFloatv(gl.GL_PROJECTION_MATRIX, self.projectionMatrix.ptr)
	self.modelViewProjectionMatrix:mul(self.projectionMatrix, self.modelViewMatrix)
end

function FrustumView:mousePan(dx, dy, screenWidth, screenHeight)
-- [[ frustum pan
	local magn = math.sqrt(dx * dx + dy * dy)
	local fdx = dx / magn
	local fdy = dy / magn
	local rotation = quatd():fromAngleAxis(-fdy, -fdx, 0, magn)
	self.angle = (self.angle * rotation):normalize()
--]]
--[[ ortho pan
	local ar = screenWidth / screenHeight
	self.pos = self.pos + vec3d(
		-dx * ar / screenWidth * self.dist,
		dy / screenHeight * self.dist,
		0)
--]]
end

function FrustumView:mouseZoom(dx, dy)
	self.dist = self.dist * math.exp(dy * -.03)
end

return FrustumView 

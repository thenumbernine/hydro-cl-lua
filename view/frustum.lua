local class = require 'ext.class'
local gl = require 'ffi.OpenGL'
local vec3 = require 'vec.vec3'
local quat = require 'vec.quat'

local FrustumView = class()

function FrustumView::init()
	self.dist = 1
	self.pos = vec3()
	self.angle = quat()
end

FrustumView.zFar = 10
FrustumView.zNear = .001
function FrustumView:projection(screenWidth, screenHeight)
	local ar = screenWidth / screenHeight
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
	gl.glRotatef(angleAxis[4], angleAxis[1], angleAxis[2], angleAxis[3])
	gl.glTranslatef(-self.pos[1], -self.pos[2], -self.pos[3])
end

function FrustumView:mousePan(dx, dy, screenWidth, screenHeight)
	local magn = math.sqrt(dx * dx + dy * dy)
	local fdx = dx / magn
	local fdy = dy / magn
	local rotation = quat():fromAngleAxis(fdy, fdx, 0, magn)
	self.angle = (rotation * self.angle):normalize()
end

function FrustumView:mouseZoom(dx, dy)
	self.dist = self.dist * math.exp(dy * -.03)
end

return FrustumView 

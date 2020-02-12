local class = require 'ext.class'
local gl = require 'gl'
local vec3d = require 'vec-ffi.vec3d'
local quatd = require 'vec-ffi.quatd'

local FrustumView = class()

function FrustumView:init()
	self.dist = 3
	self.pos = vec3d()
	self.angle = quatd(0,0,0,1)
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
	gl.glRotatef(angleAxis.w, angleAxis.x, angleAxis.y, angleAxis.z)
	gl.glTranslatef(-self.pos.x, -self.pos.y, -self.pos.z)
end

function FrustumView:mousePan(dx, dy, screenWidth, screenHeight)
	local magn = math.sqrt(dx * dx + dy * dy)
	local fdx = dx / magn
	local fdy = dy / magn
	local rotation = quatd():fromAngleAxis(fdy, fdx, 0, magn)
	self.angle = (rotation * self.angle):normalize()
end

function FrustumView:mouseZoom(dx, dy)
	self.dist = self.dist * math.exp(dy * -.03)
end

return FrustumView 

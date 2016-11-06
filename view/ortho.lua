local class = require 'ext.class'
local gl = require 'ffi.OpenGL'
local vec2 = require 'vec.vec2'

local OrthoView = class()

function OrthoView:init()
	self.zoom = vec2(.5,.5)
	self.pos = vec2()
end

function OrthoView:getOrthoBounds(ar)
	return self.pos[1] - ar * .5 / self.zoom[1], 
		self.pos[1] + ar * .5 / self.zoom[1],
		self.pos[2] - .5 / self.zoom[2],
		self.pos[2] + .5 / self.zoom[2],
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

function OrthoView:mousePan(dx, dy, screenWidth, screenHeight)
	local ar = screenWidth / screenHeight
	self.pos = self.pos + vec2(
		-dx * ar / screenWidth / self.zoom[1],
		dy / screenHeight / self.zoom[2])
end

function OrthoView:mouseZoom(dx, dy)
	self.zoom[1] = self.zoom[1] * math.exp(-dx * -.03)
	self.zoom[2] = self.zoom[2] * math.exp(dy * -.03)
end

return OrthoView 

local class = require 'ext.class'
local gl = require 'gl'
local vec2d = require 'vec-ffi.vec2d'

local OrthoView = class()

function OrthoView:init()
	self.zoom = vec2d(.5,.5)
	self.pos = vec2d()
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

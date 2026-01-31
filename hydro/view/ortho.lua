local class = require 'ext.class'
local gl = require 'gl'
local vec2d = require 'vec-ffi.vec2d'
local vec3d = require 'vec-ffi.vec3d'
local vec4x4f = require 'vec-ffi.vec4x4f'
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

	self.mvMat = vec4x4f():setIdent()
	self.projMat = vec4x4f():setIdent()
	self.mvProjMat = vec4x4f():setIdent()
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
	self.projMat:setOrtho(self:getOrthoBounds(aspectRatio))
end

function OrthoView:setupModelView()
	self.mvMat:setIdent()
end

function OrthoView:setup(aspectRatio)
	self:setupProjection(aspectRatio)
	self:setupModelView()
	self.mvProjMat:mul4x4(self.projMat, self.mvMat)
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

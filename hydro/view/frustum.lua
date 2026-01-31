local class = require 'ext.class'
local table = require 'ext.table'
local gl = require 'gl'
local vec3d = require 'vec-ffi.vec3d'
local vec4x4f = require 'vec-ffi.vec4x4f'
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

	self.mvMat = vec4x4f():setIdent()
	self.projMat = vec4x4f():setIdent()
	self.mvProjMat = vec4x4f():setIdent()
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
	self.projMat:setFrustum(self:getFrustumBounds(aspectRatio))
end

function FrustumView:setupModelView()
	local angleAxis = self.angle:toAngleAxis()
	self.mvMat
		:setTranslate(0, 0, -self.dist)
		:applyRotate(-math.rad(angleAxis.w), angleAxis.x, angleAxis.y, angleAxis.z)
		:applyTranslate(-self.pos.x, -self.pos.y, -self.pos.z)
end

function FrustumView:setup(aspectRatio)
	self:setupProjection(aspectRatio)
	self:setupModelView()
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

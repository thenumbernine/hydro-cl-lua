--[[
exact volume: 
int(theta=theta1,theta2 int(r=r1,r2 r dr) dtheta)
= 1/2 (theta2-theta1) (r2^2 - r1^2)

volume element: 
r dr dtheta
= r (r2 - r1) (theta2 - theta1)
= (r2+r1)/2 (r2 - r1) (theta2 - theta1)
= 1/2 (r2^2 - r1^2) (theta2 - theta1)

same thing, good thing
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local CoordinateSystem = require 'coord.coord'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Cylinder = class(CoordinateSystem)
Cylinder.name = 'cylinder' 

function Cylinder:init(args)
	local x, y, z = symmath.vars('x', 'y', 'z')
	self.embedded = table{x,y,z}
	
	local r, theta = symmath.vars('r', 'θ')
	self.baseCoords = table{r, theta, z}

	-- anholonomic linear transform 
	-- e_iHol = e_iHol^i partial_i 
	self.eHolToE = symmath.Matrix(
		{1, 0, 0},
		{0, 1/r, 0},
		{0, 0, 1}
	)
	
	-- https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19870019781.pdf
	-- it looks like all I need is the volume and I'm fine

	self.chart = function() return Tensor('^I', r * cos(theta), r * sin(theta), z) end
	
	Cylinder.super.init(self, args)
	
	self.rDef = (r^2 + z^2)^.5
end


function Cylinder:getCoordMapInvGLSLCode()
	local template = require 'template'	
	return template([[
vec3 coordMapInv(vec3 x) {
	//coord bounds don't seem to matter anywhere else,
	//except here in the renderer
	//TODO use the solver's bounds?
	float theta = mod(atan(x.y, x.x), <?=clnumber(2. * math.pi)?>);
	return vec3(length(x.xy), theta, x.z);
}
]], {
		clnumber = require 'cl.obj.number',
	})
end

return Cylinder

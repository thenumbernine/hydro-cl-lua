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
Cylinder.coords = {'r', 'θ', 'z'}

function Cylinder:init(args)
	local x, y, z = symmath.vars('x', 'y', 'z')
	args.embedded = table{x,y,z}
	
	local r, theta = symmath.vars('r', 'θ')

	-- [[ holonomic
	-- this requires conservation law equation delta^ij's to be swapped with g^ij's
	args.coords = table{r, theta, z}
	--]]
	--[[ anholonomic
	-- g^ij = delta^ij, so the Euler fluid equations match those in flat space
	-- however the Conn^k_jk terms still match up with the holonomic terms
	-- and the grid itself is still the holonomic coordinate system
	-- so dx's are based on the holonomic coordinate system
	-- ... other adjustments have to be made elsewhere as well 
	-- (like the extra scaling term next to the volume scale within the flux) 
	local thetaHat = symmath.var'thetaHat'
	thetaHat.base = theta
	function thetaHat:applyDiff(x) return x:diff(theta) / r end
	args.coords = table{r, thetaHat, z}
	--]]
	-- the problem with anholonomic is that 
	-- dx_at is always 1
	-- and volume is always 1
	-- so we never consider ds = r dtheta
	-- and therefore the flux never gets perturbed
	-- so TODO dx should be based on the holonomic metric?
	-- (i.e. the coordinates of the grid breakdown)?
	-- while the vector units are based on the anholonomic metri?
	-- (i.e. those with unit metric)

	-- TODO https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19870019781.pdf
	-- it looks like all I need is the volume and I'm fine

	args.chart = function() return Tensor('^I', r * cos(theta), r * sin(theta), z) end
	
	Cylinder.super.init(self, args)
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

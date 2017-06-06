local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local geometry = require 'geom.geom'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Cylinder = class(geometry)
Cylinder.name = 'cylinder' 
function Cylinder:init(args)
	local x, y, z = symmath.vars('x', 'y', 'z')
	args.embedded = table{x,y,z}:sub(1, args.solver.dim)
	
	local r, theta = symmath.vars('r', 'theta')

	--[[ holonomic
	-- this would need the conservation law equation delta^ij's swapped with g^ij's
	-- which takes a new eigen-decomposition
	args.coords = table{r, theta, z}:sub(1, args.solver.dim)
	--]]
	-- [[ anholonomic
	-- g^ij = delta^ij, so the Euler fluid equations match those in flat space
	-- however the Conn^k_jk terms still match up with the holonomic terms
	-- and the grid itself is still the holonomic coordinate system
	-- so dx's are based on the holonomic coordinate system
	local thetaHat = symmath.var'thetaHat'
	thetaHat.base = theta
	function thetaHat:applyDiff(x) return x:diff(theta) / r end
	args.coords = table{r, thetaHat, z}:sub(1, args.solver.dim)
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

	args.chart = function() 
		return ({
			function() return Tensor('^I', r) end,
			function() return Tensor('^I', r * cos(theta), r * sin(theta)) end,
			function() return Tensor('^I', r * cos(theta), r * sin(theta), z) end,
		})[args.solver.dim]()
	end
	
	Cylinder.super.init(self, args)
end

return Cylinder

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

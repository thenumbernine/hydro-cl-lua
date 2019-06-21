local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local CoordinateSystem = require 'coord.coord'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Sphere = class(CoordinateSystem)

Sphere.name = 'sphere' 
Sphere.coords = {'r', 'θ', 'φ'}

--[[
args
	volumeDim = (TODO) change volume element etc to act as if we're in a higher dimension
	
	TODO add some other arg for rearranging the coordinate order so we can do 2D simulations of θ and φ alone
--]]
function Sphere:init(args)
	args.embedded = table{symmath.vars('x', 'y', 'z')}
	local r, theta, phi = symmath.vars('r', 'θ', 'φ')

	-- [[ holonomic
	args.coords = table{r, theta, phi}
	--]]
	--[[ anholonomic
	local thetaHat = symmath.var'thetaHat'
	thetaHat.base = theta
	function thetaHat:applyDiff(expr) return expr:diff(theta) / r end

	local phiHat = symmath.var'phiHat'
	phiHat.base = phi
	function phiHat:applyDiff(expr) return expr:diff(phi) / (r * sin(theta)) end
	
	args.coords = table{thetaHat, phiHat, r}:sub(1, args.solver.dim)
	--]]
	
	args.chart = function() 
		return Tensor('^I', 
			r * sin(theta) * cos(phi), 
			r * sin(theta) * sin(phi), 
			r * cos(theta)
		) 
	end
	
	Sphere.super.init(self, args)
end

local template = require 'template'
function Sphere:getCoordMapInvGLSLCode()
	return template([[
vec3 coordMapInv(vec3 x) {
<? if solver.dim == 1 then
?>	float r = abs(x.x);
	float theta = 0.;
	float phi = 0.;
<? elseif solver.dim == 2 then	-- xy -> rθ
?>	float r = length(x.xy);
	float theta = acos(x.y / r);
	float phi = 0.;
<? elseif solver.dim == 3 then 	-- xyz - rθφ
?>	float r = length(x);
	float theta = acos(x.z / r);
	float phi = atan(x.y, x.x);
<? end 
?>	return vec3(r, theta, phi);
}
]], {
		solver = self.solver,
	})
end

return Sphere

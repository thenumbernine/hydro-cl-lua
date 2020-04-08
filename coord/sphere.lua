local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local template = require 'template'
local CoordinateSystem = require 'coord.coord'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Sphere = class(CoordinateSystem)
Sphere.name = 'sphere' 

--[[
args
	volumeDim = (TODO) change volume element etc to act as if we're in a higher dimension
	
	TODO add some other arg for rearranging the coordinate order so we can do 2D simulations of θ and φ alone
--]]
function Sphere:init(args)
	self.embedded = table{symmath.vars('x', 'y', 'z')}
	local r, theta, phi = symmath.vars('r', 'θ', 'φ')

	self.baseCoords = table{r, theta, phi}

	-- TODO same as lenExprs
	self.eHolToE = symmath.Matrix{
		{1, 0, 0},
		{0, 1/r, 0},
		{0, 0, 1/(r*symmath.sin(theta))},
	}
	
	self.chart = function() 
		return Tensor('^I', 
			r * sin(theta) * cos(phi), 
			r * sin(theta) * sin(phi), 
			r * cos(theta)
		) 
	end
	
	Sphere.super.init(self, args)

	self.vars = {
		r = r,
		x = r * sin(theta) * cos(phi),
		y = r * sin(theta) * sin(phi),
		z = r * cos(theta),
	}
end

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

function Sphere:getParallelPropagatorCode()
	return template([[

<? if coord.vectorComponent == 'holonomic' then ?>

real3 coord_parallelPropagateU0(real3 v, real3 x, real dx) {
	real rL = x.x;
	real rR = x.x + dx;
	v.y *= rL / rR;
	v.z *= rL / rR;
	return v;
}

real3 coord_parallelPropagateL0(real3 v, real3 x, real dx) {
	real rL = x.x;
	real rR = x.x + dx;
	v.y *= rR / rL;
	v.z *= rR / rL;
	return v;
}

real3 coord_parallelPropagateU1(real3 v, real3 x, real dx) {
	real r = x.x;
	real thetaL = x.y;
	real thetaR = x.y + dx;
	real sinThetaL = sin(thetaL);
	real sinThetaR = sin(thetaR);
	v.y *= r;
	v.z *= sinThetaL;
	v = real3_rotateZ(v, -dx);
	v.y /= r;
	v.z /= sinThetaR;
	return v;
}

real3 coord_parallelPropagateL1(real3 v, real3 x, real dx) {
	real r = x.x;
	real thetaL = x.y;
	real thetaR = x.y + dx;
	real sinThetaL = sin(thetaL);
	real sinThetaR = sin(thetaR);
	v.y /= r;
	v.z /= sinThetaL;
	v = real3_rotateZ(v, -dx);
	v.y *= r;
	v.z *= sinThetaR;
	return v;
}

real3 coord_parallelPropagateU2(real3 v, real3 x, real dx) {
	real r = x.x;
	real theta = x.y;
	real sinTheta = sin(theta);
	real rSinTheta = r * sinTheta;
	v.y *= r;
	v.z *= rSinTheta;
	v = real3_rotateZ(v, theta);
	v = real3_rotateX(v, -dx);
	v = real3_rotateZ(v, -theta);
	v.y /= r;
	v.z /= rSinTheta;
	return v;
}

real3 coord_parallelPropagateL2(real3 v, real3 x, real dx) {
	real r = x.x;
	real theta = x.y;
	real sinTheta = sin(theta);
	real rSinTheta = r * sinTheta;
	v.y /= r;
	v.z /= rSinTheta;
	v = real3_rotateZ(v, theta);
	v = real3_rotateX(v, -dx);
	v = real3_rotateZ(v, -theta);
	v.y *= r;
	v.z *= rSinTheta;
	return v;
}

<? else ?>

#define coord_parallelPropagateU0(v,x,dx) (v)

real3 coord_parallelPropagateU1(real3 v, real3 x, real dx) {
	return real3_rotateZ(v, -dx);
}

real3 coord_parallelPropagateU2(real3 v, real3 x, real dx) {
	return real3_rotateX(-dx);
}

<? end ?>

]], {
		coord = self,
	})
end

return Sphere

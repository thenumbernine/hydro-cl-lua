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
local template = require 'template'	
local CoordinateSystem = require 'hydro.coord.coord'

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

	self.vars = {
		r = (r^2 + z^2)^.5,
		x = r * cos(theta),
		y = r * sin(theta),
		z = z,
	}
end


function Cylinder:getCoordMapInvGLSLCode()
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

function Cylinder:getParallelPropagatorCode()
	return template([[

#define coord_parallelPropagateU2(v,x,dx) (v)

<? if coord.vectorComponent == 'holonomic' then ?>

#define coord_parallelPropagateL2 coord_parallelPropagateU2

real3 coord_parallelPropagateU0(real3 v, real3 x, real dx) {
	real rL = x.x;
	real rR = x.x + dx;
	v.y *= rL / rR;
	return v;
}

real3 coord_parallelPropagateL0(real3 v, real3 x, real dx) {
	real rL = x.x;
	real rR = x.x + dx;
	v.y *= rR / rL;
	return v;
}

real3 coord_parallelPropagateU1(real3 v, real3 x, real dx) {
	real rL = x.x;
	real rR = x.x + dx;
	v.y *= rL;
	v = real3_rotateZ(v, -dx);
	v.y /= rR;
	return v;
}

real3 coord_parallelPropagateL1(real3 v, real3 x, real dx) {
	real rL = x.x;
	real rR = x.x + dx;
	v.y /= rL;
	v = real3_rotateZ(v, -dx);
	v.y *= rR;
	return v;
}

<? elseif coord.vectorComponent == 'anholonomic' then ?>

#define coord_parallelPropagateU0(v,x,dx) (v)

// propagate a vector's components
// anholonomic, meaning no rescaling
real3 coord_parallelPropagateU1(real3 v, real3 x, real dx) {
	return real3_rotateZ(v, -dx);
}

<? end ?>
]], {
		coord = self,
	})
end

return Cylinder

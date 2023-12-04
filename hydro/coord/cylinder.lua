--[[
exact volume:
int(φ=φ1,φ2 int(r=r1,r2 r dr) dφ)
= 1/2 (φ2-φ1) (r2^2 - r1^2)

volume element:
r dr dφ
= r (r2 - r1) (φ2 - φ1)
= (r2+r1)/2 (r2 - r1) (φ2 - φ1)
= 1/2 (r2^2 - r1^2) (φ2 - φ1)

same thing, good thing
--]]

local table = require 'ext.table'
local symmath = require 'symmath'
local template = require 'template'
local CoordinateSystem = require 'hydro.coord.coord'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Cylinder = CoordinateSystem:subclass()
Cylinder.name = 'cylinder'

function Cylinder:init(args)
	local x, y, z = symmath.vars('x', 'y', 'z')
	self.embedded = table{x,y,z}

	-- this r is the 2D r, not the 3D r ... I should use a different letter ...
	-- Wiki uses 'rho' to separate from 3D r
	-- but my spherical sinh-remapped coordinate system also uses rho ...
	-- ... should I change that one to 'eta' and change this to 'rho'?
	-- ... or should I just not bother to keep my local Coord variable names unique?
	local r2, phi = symmath.vars('r', 'φ')
	self.baseCoords = table{r2, phi, z}

	-- anholonomic linear transform
	-- e_iHol = e_iHol^i partial_i
	self.eHolToE = symmath.Matrix(
		{1, 0, 0},
		{0, 1/r2, 0},
		{0, 0, 1}
	)

	-- https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19870019781.pdf
	-- it looks like all I need is the volume and I'm fine

	self.chart = function()
		return Tensor('^I',
			r2 * cos(phi),
			r2 * sin(phi),
			z
		)
	end

	Cylinder.super.init(self, args)

	self.vars = {
		r = (r2^2 + z^2)^.5,
		theta = symmath.atan2(r2, z),
		phi = phi,
		x = r2 * cos(phi),
		y = r2 * sin(phi),
		z = z,
	}
end


function Cylinder:getCoordMapInvModuleCode()
	return [[
real3 coordMapInv(real3 x) {
	//coord bounds don't seem to matter anywhere else,
	//except here in the renderer
	//TODO use the solver's bounds?
	real phi = fmod(fmod(atan2(x.y, x.x), 2. * M_PI) + 2. * M_PI, 2. * M_PI);
	return real3(sqrt(x.x*x.x + x.y*x.y), phi, x.z);
}
]]
end

function Cylinder:getModuleDepends_coord_parallelPropagate()
	return {'rotate'}
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

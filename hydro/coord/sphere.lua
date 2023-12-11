local table = require 'ext.table'
local symmath = require 'symmath'
local template = require 'template'
local CoordinateSystem = require 'hydro.coord.coord'

local sin, cos = symmath.sin, symmath.cos
local Tensor = symmath.Tensor

local Sphere = CoordinateSystem:subclass()
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
	self.eHolToE = symmath.Matrix(
		{1, 0, 0},
		{0, 1/r, 0},
		{0, 0, 1/(r*symmath.sin(theta))}
	)

	self.chart = function()
		return Tensor('^I',
			r * sin(theta) * cos(phi),
			r * sin(theta) * sin(phi),
			r * cos(theta)
		)
	end

	Sphere.super.init(self, args)

	self.vars = {
		x = r * sin(theta) * cos(phi),
		y = r * sin(theta) * sin(phi),
		z = r * cos(theta),
		r = r,
		theta = theta,
		phi = phi,
	}
end

function Sphere:getCoordMapInvModuleCode()
	return template([[
real3 coordMapInv(real3 x) {
<? if solver.dim == 1 then
?>	real r = fabs(x.x);
	real theta = 0.;
	real phi = 0.;
<? elseif solver.dim == 2 then	-- xy -> rθ
?>	real r = sqrt(x.x*x.x + x.y*x.y);
	real theta = acos(x.y / r);
	real phi = 0.;
<? elseif solver.dim == 3 then 	-- xyz - rθφ
?>	real r = length(x);
	real theta = acos(x.z / r);
	real phi = atan2(x.y, x.x);
<? end
?>	return _real3(r, theta, phi);
}
]], {
		solver = self.solver,
	})
end

function Sphere:getModuleDepends_coord_parallelPropagate()
	return {'rotate'}
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
	return real3_rotateX(v, -dx);
}

<? end ?>

]], {
		coord = self,
	})
end

function Sphere:createCellStruct()
	Sphere.super.createCellStruct(self)
	self.cellStructFields:insert{name='r', type='real'}
end

function Sphere:fillGridCellBuf(cellsCPU)
	local solver = self.solver

	local r, theta, phi = self.baseCoords:unpack()
	local calcR, code = symmath.export.Lua:toFunc{
		output = {self.vars.r},
		input = {{r=r}, {theta=theta}, {phi=phi}},
	}

	local index = 0
	for k=0,tonumber(solver.gridSize.z)-1 do
		local phi = solver.dim >= 3
			and ((k + .5 - solver.numGhost) / (tonumber(solver.gridSize.z) - 2 * solver.numGhost) * (solver.maxs.z - solver.mins.z) + solver.mins.z)
			or (.5 * (solver.maxs.z + solver.mins.z))
		for j=0,tonumber(solver.gridSize.y)-1 do
			local theta = solver.dim >= 2
				and ((j + .5 - solver.numGhost) / (tonumber(solver.gridSize.y) - 2 * solver.numGhost) * (solver.maxs.y - solver.mins.y) + solver.mins.y)
				or (.5 * (solver.maxs.y + solver.mins.y))
			for i=0,tonumber(solver.gridSize.x)-1 do
				local r = solver.dim >= 1
					and ((i + .5 - solver.numGhost) / (tonumber(solver.gridSize.x) - 2 * solver.numGhost) * (solver.maxs.x - solver.mins.x) + solver.mins.x)
					or (.5 * (solver.maxs.x + solver.mins.x))
				cellsCPU[index].pos.x = r
				cellsCPU[index].pos.y = theta
				cellsCPU[index].pos.z = phi
				cellsCPU[index].r = calcR(r, theta, phi)
				index = index + 1
			end
		end
	end
end

return Sphere

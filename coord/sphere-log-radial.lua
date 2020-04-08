local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local CoordinateSystem = require 'coord.coord'

local sin, cos = symmath.sin, symmath.cos
local sinh = symmath.sinh
local Tensor = symmath.Tensor

local SphereLogRadial = class(CoordinateSystem)
SphereLogRadial.name = 'sphere-log-radial'

--[[
args
	volumeDim = (TODO) change volume element etc to act as if we're in a higher dimension
	
	TODO add some other arg for rearranging the coordinate order so we can do 2D simulations of θ and φ alone
--]]
function SphereLogRadial:init(args)
	self.embedded = table{symmath.vars('x', 'y', 'z')}
	local rho, theta, phi = symmath.vars('ρ', 'θ', 'φ')

	self.baseCoords = table{rho, theta, phi}

	-- 2017 Ruchlin, Etienne after eqn 42
	local solver = args.solver
	local amplitude = 1000
	local sinh_w = .15
self.amplitude = amplitude
self.sinh_w = sinh_w
	local rmax = solver.maxs.x
	local r = symmath.var('r', {rho})

	local rDef = amplitude * sinh(rho / sinh_w) / math.sinh(1 / sinh_w)

	self.vars = {
		x = rDef * sin(theta) * cos(phi),
		y = rDef * sin(theta) * sin(phi),
		z = rDef * cos(theta),
		r = rDef,
		theta = theta,
		phi = phi,
	}

	local r_for_rho = rDef
	self.replvars = table{
		{r:diff(rho, rho, rho), r_for_rho:diff(rho, rho, rho)()},
		{r:diff(rho, rho), r_for_rho:diff(rho, rho)()},
		{r:diff(rho), r_for_rho:diff(rho)()},
		{r, r_for_rho},
	}

	if cmdline.coordVerbose then
		print(r:eq(r_for_rho))
	end
	self.rho_for_r = r:eq(r_for_rho):solve(rho):rhs()
	if cmdline.coordVerbose then
		print(rho:eq(self.rho_for_r))
	end
	self.eHolToE = symmath.Matrix{
		{1/r:diff(rho), 0, 0},
		{0, 1/r, 0},
		{0, 0, 1/(r*symmath.sin(theta))},
	}
	
	self.chart = function() 
		return Tensor('^I', 
			r * symmath.sin(theta) * symmath.cos(phi), 
			r * symmath.sin(theta) * symmath.sin(phi), 
			r * symmath.cos(theta)
		) 
	end
	
	SphereLogRadial.super.init(self, args)
end

local template = require 'template'
function SphereLogRadial:getCoordMapInvGLSLCode()
	return template([[

float sinh(float x) { return .5 * (exp(x) - exp(-x)); }
float cosh(float x) { return .5 * (exp(x) + exp(-x)); }
float asinh(float x) { return log(x + sqrt(x*x + 1.)); }

vec3 coordMapInv(vec3 pt) {
<? if solver.dim == 1 then
?>	float r = abs(pt.x);
	float theta = 0.;
	float phi = 0.;
<? elseif solver.dim == 2 then	-- xy -> rθ
?>	float r = length(pt.xy);
	float theta = acos(pt.y / r);
	float phi = 0.;
<? elseif solver.dim == 3 then 	-- xyz - rθφ
?>	float r = length(pt);
	float theta = acos(pt.z / r);
	float phi = atan(pt.y, pt.x);
<? end 
?>	
	float rho = <?=coord:compile(coord.rho_for_r)?>;
	return vec3(rho, theta, phi);
}
]], {
		symmath = require 'symmath',
		coord = self,
		solver = self.solver,
	})
end

function SphereLogRadial:getParallelPropagatorCode()
	return template([[
<? local clnumber = require 'cl.obj.number' ?>

<? if coord.vectorComponent == 'holonomic' then ?>

real3 coord_parallelPropagateU0(real3 v, real3 x, real dx) {
	const real w = <?=clnumber(coord.sinh_w)?>;
	real rhoL = x.x;
	real rhoR = x.x + dx;
	real coshLR = cosh(rhoL/w) / cosh(rhoR/w);
	real sinhLR = sinh(rhoL/w) / sinh(rhoR/w);
	v.x *= coshLR;
	v.y *= sinhLR;
	v.z *= sinhLR;
	return v;
}

real3 coord_parallelPropagateL0(real3 v, real3 x, real dx) {
	const real w = <?=clnumber(coord.sinh_w)?>;
	real rhoL = x.x;
	real rhoR = x.x + dx;
	real coshRL = cosh(rhoR/w) / cosh(rhoL/w);
	real sinhRL = sinh(rhoR/w) / sinh(rhoL/w);
	v.x *= coshRL;
	v.y *= sinhRL;
	v.z *= sinhRL;
	return v;
}

real3 coord_parallelPropagateU1(real3 v, real3 x, real dx) {
	const real w = <?=clnumber(coord.sinh_w)?>;
	real rho = v.x;
	real thetaL = x.y;
	real thetaR = x.y + dx;
	real s = w * sinh(rho/w) / cosh(rho/w);
	v.y *= s;
	v = real3_rotateZ(v, -dx);
	v.y /= s;
	v.z *= fabs(sin(thetaL) / sin(thetaR));
	return v;
}

real3 coord_parallelPropagateL1(real3 v, real3 x, real dx) {
	const real w = <?=clnumber(coord.sinh_w)?>;
	real rho = v.x;
	real thetaL = x.y;
	real thetaR = x.y + dx;
	real s = w * sinh(rho/w) / cosh(rho/w);
	v.y /= s;
	v = real3_rotateZ(v, -dx);
	v.y *= s;
	v.z *= fabs(sin(thetaR) / sin(thetaL));
	return v;
}

// TODO here ... fix these
real3 coord_parallelPropagateU2(real3 v, real3 x, real dx) {
	const real w = <?=clnumber(coord.sinh_w)?>;
	real rho = x.x;
	real s = w * sinh(rho/w) / cosh(rho/w);
	real theta = x.y;
	real sinTheta = sin(theta);
	real sSinTheta = s * sinTheta;
	v.y *= s;
	v.z *= sSinTheta;
	v = real3_rotateZ(v, theta);
	v = real3_rotateX(v, -dx);
	v = real3_rotateZ(v, -theta);
	v.y /= s;
	v.z /= sSinTheta;
	return v;
}

real3 coord_parallelPropagateL2(real3 v, real3 x, real dx) {
	const real w = <?=clnumber(coord.sinh_w)?>;
	real rho = x.x;
	real s = w * sinh(rho/w) / cosh(rho/w);
	real theta = x.y;
	real sinTheta = sin(theta);
	real sSinTheta = s * sinTheta;
	v.y /= s;
	v.z /= sSinTheta;
	v = real3_rotateZ(v, theta);
	v = real3_rotateX(v, -dx);
	v = real3_rotateZ(v, -theta);
	v.y *= s;
	v.z *= sSinTheta;
	return v;
}

<? else 
	error "still need to do anholonomic"
end ?>

]], {
		coord = self,
	})
end

return SphereLogRadial

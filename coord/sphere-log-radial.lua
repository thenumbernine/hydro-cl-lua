local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local CoordinateSystem = require 'coord.coord'

local sin, cos = symmath.sin, symmath.cos
local sinh = symmath.sinh
local Tensor = symmath.Tensor

local Sphere = class(CoordinateSystem)
Sphere.name = 'sphere-log-radial'

--[[
args
	volumeDim = (TODO) change volume element etc to act as if we're in a higher dimension
	
	TODO add some other arg for rearranging the coordinate order so we can do 2D simulations of θ and φ alone
--]]
function Sphere:init(args)
	self.embedded = table{symmath.vars('x', 'y', 'z')}
	local rho, theta, phi = symmath.vars('ρ', 'θ', 'φ')

	self.baseCoords = table{rho, theta, phi}

	-- 2018 Ruchlin, Etienne
	local solver = args.solver
	local amplitude = 1000
	local sinh_w = .15
self.amplitude = amplitude
self.sinh_w = sinh_w
	local rmax = solver.maxs.x
	local r = symmath.var('r', {rho})

	local rDef = amplitude * sinh(rho / sinh_w) / math.sinh(1 / sinh_w)
self.rDef = rDef

	self.vars = {
		r = rDef,
		x = rDef * sin(theta) * cos(phi),
		y = rDef * sin(theta) * sin(phi),
		z = rDef * cos(theta),
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
	
	Sphere.super.init(self, args)
end

local template = require 'template'
function Sphere:getCoordMapInvGLSLCode()
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

return Sphere

local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local CoordinateSystem = require 'hydro.coord.coord'
local template = require 'template'

local var = symmath.var
local frac = symmath.frac
local sin, cos = symmath.sin, symmath.cos
local sinh = symmath.sinh
local Tensor = symmath.Tensor

local SphereLogRadial = class(CoordinateSystem)
SphereLogRadial.name = 'sphere-sinh-radial'

-- I was trying to hold off on this, thinking it would save something somewhere, because the subsequent conn calcs were slow
-- but I don't need all them anyways (cue code module system)
-- and it looks like it was stiff-arming me when it came to bssnok-fd-sym analytic simplifications.
-- if we don't defer then calculating connections is incredibly slow
-- bssnok-fd-sym needs it not defered.  TODO substitute it in bssnok-fd-sym

--[[
args
	volumeDim = (TODO) change volume element etc to act as if we're in a higher dimension
	
	TODO add some other arg for rearranging the coordinate order so we can do 2D simulations of θ and φ alone
	amplitude
	sinh_w
--]]
function SphereLogRadial:init(args)
	self.embedded = table{symmath.vars('x', 'y', 'z')}
	local rho, theta, phi = symmath.vars('ρ', 'θ', 'φ')

	self.baseCoords = table{rho, theta, phi}

	-- 2017 Ruchlin, Etienne after eqn 42
	local solver = args.solver
	self.amplitude = args.amplitude or 1000 
	self.sinh_w = args.sinh_w or .15
--	local rmax = solver.maxs.x	-- not used (also not initialized at this point)

	local amplitude, sinh_w

	-- and now that all the coord code now depends on the compile-time solver vars, lets make sure everyone can see them
	-- (technically anyone who uses any code that comes from solver.coord)
	-- solver.sharedModulesEnabled['eqn.guiVars.compileTime'] = true
	-- but this isn't defined yet in solver, so move it later ...
	self.amplitude_var = var'AMPL'
	self.sinh_w_var = var'SINHW'
	amplitude = self.amplitude_var
	sinh_w = self.sinh_w_var
	
	-- TODO repl vars for solver->coord_sinh_w and solver->coord_amplitude, and make these coord parameters?
	-- or they can just be compile-time, but both can be accomplished with solver's guiVars
	
	local rDef = amplitude * sinh(rho / sinh_w) / sinh(frac(1, sinh_w))

	local r = symmath.var('r', {rho})
	
	self.vars = {
		x = rDef * sin(theta) * cos(phi),
		y = rDef * sin(theta) * sin(phi),
		z = rDef * cos(theta),
		r = rDef,
		theta = theta,
		phi = phi,
	}

	-- don't replace these until we are compiling or integrating, to save on simplification time
	local r_for_rho = rDef
	self.repls = table{
		r:diff(rho, rho, rho):eq(r_for_rho:diff(rho, rho, rho)()),
		r:diff(rho, rho):eq(r_for_rho:diff(rho, rho)()),
		r:diff(rho):eq(r_for_rho:diff(rho)()),
		r:eq(r_for_rho),
	}

	if cmdline.coordVerbose then
		print(r:eq(r_for_rho))
	end
	self.rho_for_r = r:eq(r_for_rho):solve(rho):rhs()
	if cmdline.coordVerbose then
		print(rho:eq(self.rho_for_r))
	end

	self.eHolToE = symmath.Matrix(
		{1/r:diff(rho), 0, 0},
		{0, 1/r, 0},
		{0, 0, 1/(r*symmath.sin(theta))}
	)
	
	self.chart = function() 
		return Tensor('^I', 
			r * symmath.sin(theta) * symmath.cos(phi), 
			r * symmath.sin(theta) * symmath.sin(phi), 
			r * symmath.cos(theta)
		) 
	end
	
	SphereLogRadial.super.init(self, args)
end

function SphereLogRadial:getModuleDepends_coordMap() 
	return {'eqn.guiVars.compileTime'}
end

-- TODO only for GLSL, depend on 'sinh'
-- I could put '#ifdef GLSL' around the sinh and cosh modules, but what about asinh which C needs as well?
-- I could put a 'glsl_code' option for each module, and give them to the sinh,cosh,acosh, but ... same question?
-- I could add them as glslModulesEnabled[] ... but what would guarantee the order of code generation?
-- I could add/remove sinh from coordMap's depends or code before/after adding it to the GLSL ... kind of ugly and very specific to sphere-sinh-radial only
-- ... see this is a good argument for replacing the 'depends' system with a code gen execution and an 'include' function with 'pragma once' ability
--		then coordMap for sphere-sinh-radial could have code that said "if included from glsl then include sinh end"
--		also you can put the compile() symbolic codegen into the same code doing the include's of dependencies
-- 		so to do this, get rid of 'depends', replace 'code', 'header', 'typecode' strings with functions, and give them an 'include' function that pulls in from other modules by their name
-- for now I'll just make a separate 'coordMapGLSL' just like 'coordMap' but with the GLSL depends
function SphereLogRadial:getModuleDepends_coordMapGLSL() 
	return {
		'sinh',
		-- This is in sharedModules ... so should sharedModules be included in all GLSL modules?
		-- until I decide, I'll just put this entry for GLSL here:
		'eqn.guiVars.compileTime',	-- for AMPL and SINHW
	}
end
function SphereLogRadial:getModuleDepends_coordMapInvGLSL() 
	return {
		'sinh',
		'cosh',
		'asinh',
		'eqn.guiVars.compileTime',	-- for AMPL and SINHW
	}
end

function SphereLogRadial:initCodeModules(...)
	SphereLogRadial.super.initCodeModules(self, ...)
	-- since now all coord code uses AMPL and SINHW macros ...
	-- in fact to have it first ...
	-- which means TODO it should be in the code module 'depends' of all code that 'coord' generates
	-- so how about a 'getDepends' that is subclassed?
	self.solver.eqn:addGuiVar{name='AMPL', value=self.amplitude, compileTime=true}
	self.solver.eqn:addGuiVar{name='SINHW', value=self.sinh_w, compileTime=true}
	self.solver.sharedModulesEnabled['eqn.guiVars.compileTime'] = true
end

--[[
This function is especially used in GLSL, and GLSL doesn't like converting literal integers to floats, so use '.'s after all your numbers
--]]
function SphereLogRadial:getCoordMapInvModuleCode()
	return template([[
real3 coordMapInv(real3 pt) {
	//this part matches sphere ... hmm
<? if solver.dim == 1 then
?>	real r = fabs(pt.x);
	real theta = .5*M_PI;
	real phi = 0.;
<? elseif solver.dim == 2 then	-- xy -> rθ
?>	real r = real3_len(pt.xy);
	real theta = acos(pt.y / r);
	real phi = 0.;
<? elseif solver.dim == 3 then 	-- xyz - rθφ
?>	real r = real3_len(pt);
	real theta = acos(pt.z / r);
	real phi = atan2(pt.y, pt.x);
	if (phi < 0.) phi += 2. * M_PI;
<? end 
?>	

	real rho = <?=
-- the default expression uses 'pt' for the input arg, and pt.x for 'r' ...
-- ... but we're already using 'pt' for xyz and 'pt.r' for x ...
-- and if I replace it with var'r' ...
-- ... then will coord:compile use coord.repls to replace that with r_for_rho?
--		coord:compile(coord.rho_for_r:replace(coord.baseCoords[1], require 'symmath'.var'r'))
		require 'symmath.export.C'(coord.rho_for_r
			:replace(coord.amplitude_var, coord.amplitude)
			:replace(coord.sinh_w_var, coord.sinh_w)
			:replace(coord.baseCoords[1], require 'symmath'.var'r')
		)
	?>;
	
	if (rho == 0. || theta == 0. || theta == M_PI) return _real3(0.,0.,0.);
	return _real3(rho, theta, phi);
}
]], {
		symmath = require 'symmath',
		coord = self,
		solver = self.solver,
	})
end

function SphereLogRadial:getModuleDepends_coord_parallelPropagate()
	return {
		'rotate',
	}
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

<? else ?>

#define coord_parallelPropagateU0(v, x, dx) (v)

real3 coord_parallelPropagateU1(real3 v, real3 x, real dx) {
	return real3_rotateZ(v, -dx);
}

real3 coord_parallelPropagateU2(real3 v, real3 x, real dx) {
	return real3_rotateX(v, -dx);
}

<? end ?>

]], {
		coord = self,
		solver = self.solver,
	})
end

function SphereLogRadial:createCellStruct()
	SphereLogRadial.super.createCellStruct(self)
	self.cellStruct.vars:insert{name='r', type='real'}
end

function SphereLogRadial:fillGridCellBuf(cellsCPU)
	local solver = self.solver

	local symmath = require 'symmath'
	local rho, theta, phi = self.baseCoords:unpack()
	local calcR, code = symmath.export.Lua:toFunc{
		output = {
			self.vars.r
				:replace(self.amplitude_var, self.amplitude)
				:replace(self.sinh_w_var, self.sinh_w)
		},
		input = {{rho=rho}, {theta=theta}, {phi=phi}},
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
				local rho = solver.dim >= 1
					and ((i + .5 - solver.numGhost) / (tonumber(solver.gridSize.x) - 2 * solver.numGhost) * (solver.maxs.x - solver.mins.x) + solver.mins.x)
					or (.5 * (solver.maxs.x + solver.mins.x))
				cellsCPU[index].pos.x = rho
				cellsCPU[index].pos.y = theta
				cellsCPU[index].pos.z = phi
				cellsCPU[index].r = calcR(rho, theta, phi)
				index = index + 1
			end
		end
	end
end

return SphereLogRadial

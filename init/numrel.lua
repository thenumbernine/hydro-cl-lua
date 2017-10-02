local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local template = require 'template'
local vec3 = require 'vec.vec3'
local clnumber = require 'cl.obj.number'
local InitCond = require 'init.init'

require 'common'(_G)

local function getTemplateEnv(solver)
	return {
		solver = solver,
		xNames = xNames,
		symNames = symNames,
		from3x3to6 = from3x3to6,
		from6to3x3 = from6to3x3,
		clnumber = clnumber,
	}
end

local function compileC(expr, name, vars)
	assert(type(expr) == 'table', "expected table, found "..type(expr))
	if symmath.Expression.is(expr) then 
		expr = expr()
	
		local range = require 'ext.range'
		-- replace pow(x,2) with x*x
		expr = expr:map(function(x)
			if symmath.op.pow.is(x) 
			and symmath.Constant.is(x[2])
			then
				local value = assert(x[2].value)
				if value > 0 and value == math.floor(value) then
					if value == 1 then
						return x[1]
					else
						return symmath.op.mul(range(value):map(function() 
							return symmath.clone(x[1])
						end):unpack())
					end
				end
			end
		end)		
		
		print('compiling '..name..':')
		print(expr)
		local code = expr:compile(vars, 'C'), name
	
		-- ugh...
		code = code:gsub('sqrt%(', '(real)sqrt((real)')
		
		print(code)
		return code
	end
	return table.map(expr, function(v,k) return compileC(v,k,vars) end), name
end

local NumRelInitCond = class(InitCond)

function NumRelInitCond:refreshInitStateProgram(solver)
	NumRelInitCond.super.refreshInitStateProgram(self, solver)
	solver.initDerivsKernelObj = solver.initStateProgramObj:kernel('initDerivs', solver.UBuf)
end

function NumRelInitCond:resetState(solver)
	NumRelInitCond.super.resetState(self, solver)
	solver:boundary()
	solver.initDerivsKernelObj()
end

function NumRelInitCond:getCodePrefix(solver)
	-- looks like all num rel solvers might need this
	-- maybe I should put it in InitCond?
	
	local alphaVar = symmath.var'alpha'
	-- TODO each eqn must be stated as a guiVar of the solver
	-- make a parent class or something for all num rel eqns
	local fGuiVar = solver.eqn.guiVars.f
	local fLuaCode = fGuiVar.options[fGuiVar.value]
	
	local f = assert(loadstring([[
local alpha, symmath = ...
local log = symmath.log
return ]]..fLuaCode))(alphaVar, symmath)
	f = symmath.clone(f)
	local dalpha_f = f:diff(alphaVar)()

	local codes = table()
	codes.f = compileC(f, 'f', {alphaVar})
	codes.dalpha_f = compileC(dalpha_f, 'dalpha_f', {alphaVar})

	return codes:map(function(code,name,t)
		return 'real calc_'..name..code, #t+1
	end):concat'\n'		
end

-- I'm working on a unified initial condition code for 1D and 3D NR problems:
--[[
args:
	solver = the solver
	getCodes = function that accepts a table of the expressions alpha, beta, gamma, K
				and returns which expressions are needed by this particular equation

	args = {x,y,z} spatial basis variables
	alpha = lapse expression
	beta = shift expressions 
	gamma = 3-metric
	K = extrinsic curvature

	symmetric 3x3 matrices (particularly gamma & K) are stored as {xx,xy,xz,yy,yz,zz}

	(still working on)
	density, pressure = lua function

	preCompile = function to call on expressions before compile
--]]
local function initNumRel(args)	
	local vars = assert(args.vars)
	
	local exprs = table{
		alpha = assert(args.alpha),
		beta = {table.unpack(args.beta or {0,0,0})},	-- optional
		gamma = {table.unpack(args.gamma)},
		K = {table.unpack(args.K)}
	}
	assert(#exprs.beta == 3, "didn't find 3 beta vars")
	assert(#exprs.gamma == 6, "didn't find 6 gamma vars")
	assert(#exprs.K == 6, "didn't find 6 K vars")

	local function toExpr(expr, name)
		if type(expr) == 'number' then expr = symmath.Constant(expr) end
		if type(expr) == 'table' then
			if not expr.isa then
				expr = table.map(expr, toExpr)
			end
			-- simplify?
			if symmath.Expression.is(expr) then
				expr = expr()
			end
		end
		return expr, name
	end
	print('converting everything to expressions...')
	exprs = table.map(exprs, toExpr)
	print('...done converting everything to expressions')
	
	-- ADM
	
	-- assuming this is one of the ADM 1D solvers
	print('compiling expressions...')
	local codes = table.map(
		assert(args.getCodes(exprs, vars, args), "getCodes needs to return something"),
		function(v,k) 
			if args.preCompile then
				v = args.preCompile(v)
			end
			return compileC(v,k,vars) 
		end)
	print('...done compiling expressions')
	
	-- here's the lapse conditions 

	local alphaVar = symmath.var'alpha'
	-- TODO each eqn must be stated as a guiVar of the solver
	-- make a parent class or something for all num rel eqns
	local fGuiVar = args.solver.eqn.guiVars.f
	local fLuaCode = fGuiVar.options[fGuiVar.value]
	
	local f = assert(loadstring([[
local alpha, symmath = ...
local log = symmath.log
return ]]..fLuaCode))(alphaVar, symmath)
	f = symmath.clone(f)
	local dalpha_f = f:diff(alphaVar)()

	codes.f = compileC(f, 'f', {alphaVar})
	codes.dalpha_f = compileC(dalpha_f, 'dalpha_f', {alphaVar})

	return table.map(codes, function(code,name,t)
		return 'real calc_'..name..code, #t+1
	end):concat'\n'
end

function NumRelInitCond:buildFCCode(solver, diff)
	local alphaVar = symmath.var'alpha'
	local fGuiVar = solver.eqn.guiVars.f
	local fLuaCode = fGuiVar.options[fGuiVar.value]
	
	local f = assert(loadstring([[
local alpha, symmath = ...
local log = symmath.log
return ]]..fLuaCode))(alphaVar, symmath)
	f = symmath.clone(f)
	local fCCode = compileC(f, 'f', {alphaVar})

	return fCCode
end

return table{
	-- from 1997 Alcubierre "The appearance of coorindate shocks in hyperbolic formalisms of General Relativity".
	{
		name = 'gaussian perturbation',
		init = function(self, solver)
			local size = solver.maxs[1] - solver.mins[1]
			
			-- 1997 Alcubierre uses amplitude of 5 on a grid size of 300 with dx=1
			--local H = 5 / 300 * size
			--local sigma = 10 / 300 * size
			-- 1998 Bona et al use amplitude of .05 on a grid size of 100 with dx=.01
			--local H = .05 * size
			--local sigma = math.sqrt(.05) * size
			-- messing with it on my own
			local H = .01 * size
			local sigma = .1 * size

			-- TODO if only OpenCL allowed something like uniforms ...
			solver.eqn:addGuiVars{
				{name = 'H', value = H},
				{name = 'sigma', value = sigma},
			}
		end,
		initState = function(self, solver)
			return [[
	
	real3 c = real3_sub(x, mids);
	real s = real3_lenSq(c);

	const real H = gui_H;
	const real sigma = gui_sigma;
	const real sigma2 = sigma * sigma;
	const real sigma4 = sigma2 * sigma2;
	real h = H * exp(-s / sigma2);

	//h,i = (H exp(-(x-c)^2 / sigma^2)),i
	// = -2 h / sigma^2 (x_i-c_i)
	real3 dh = real3_scale(real3_sub(x, mids), -2. * h / sigma2);

	sym3 delta_ll = sym3_ident();

	//h,ij = h,i,j = (-2 h / sigma^2 (x_i-c_i)),j
	// = -2 h,j / sigma^2 (x_i-c_i) - 2 h delta_ij / sigma^2
	// = -2 (-2 h / sigma^2 (x-c)_i (x-c)_j) / sigma^2 - 2 h delta_ij / sigma^2
	// = 4 h (x-c)_i (x-c)_j / sigma^4 - 2 h delta_ij / sigma^2
	sym3 d2h = sym3_sub(
		sym3_scale(real3_outer(c, c), 4. * h / sigma4),
		sym3_scale(sym3_ident(), 2. * h / sigma2));

#if 0
	alpha = .5 * (
		(1 + sqrt(1 + kappa)) * sqrt((1-dh.x)/(1+dh.x))
		- kappa / (1 + sqrt(1 + kappa)) * sqrt((1+dh.x)/(1-dh.x))
	);
#endif

	gamma_ll = sym3_sub(delta_ll, real3_outer(dh, dh));
	K_ll = sym3_scale(d2h, -1./sqrt(1. - real3_lenSq(dh)));

//enable this if you want the ADM 3D run in 1D to match the ADM 1D's
//disable this if you want things to run in higher dimensions
#if dim == 1
	gamma_ll = _sym3(gamma_ll.xx, 0,0,1,0,1);
	K_ll = _sym3(K_ll.xx, 0,0,0,0,0);
#endif
]]
		end,
	},
	-- from 2012 Alic et. al. "Conformal and covariant formulations of the Z4 system with constraint-violation damping"
	{
		name = 'plane gauge wave',
		init = function(self, solver)
			solver.eqn:addGuiVars{
				{name = 'A', value = .1},
				{name = 'L', value = 1},
			}
		end,
		initState = function(self, solver)
			return [[
	real h = 1. - gui_A * sin((2. * M_PI / gui_L) * x.x);
	alpha = sqrt(h);
	gamma_ll.xx = h;

]]
		end,
	},
	{
		name = 'Alcubierre warp bubble',
		init = function(self, solver)
			solver.eqn:addGuiVars{
				{name = 'R', value = .5},		-- warp bubble radius
				{name = 'sigma', value = 8},	-- warp bubble thickness
				{name = 'speed', value = .1},	-- warp bubble speed
			}
		end,
		initState = function(self, solver)
			return [[
	real x_s = 0;	//speed * t
	real v_s = gui_speed;
	
	real3 y = x; y.x -= x_s;
	real r_s = real3_len(y);

#define cosh(x)		(.5 * (exp(x) + exp(-x)))
#define dtanh(x) 	(1./(cosh(x)*cosh(x)))

	real fnum = tanh(gui_sigma * (r_s + gui_R)) - tanh(gui_sigma * (r_s - gui_R));
	real fdenom = 2 * tanh(gui_sigma * gui_R);
	real f = fnum / fdenom;

	beta_u.x = -v_s * f;

	real3 dx_r_s = real3_scale(y, 1. / r_s);

	real3 dx_f = real3_scale(dx_r_s, 
		gui_sigma * (
			dtanh(gui_sigma * (r_s + gui_R)) 
			- dtanh(gui_sigma * (r_s - gui_R))
		) / fdenom);

	alpha = 1;

	K_ll.xx = -v_s * dx_f.x / alpha;
	K_ll.xy = -v_s * dx_f.y / (2 * alpha);
	K_ll.xz = -v_s * dx_f.z / (2 * alpha);
]]
		end,
	},
	{	-- take the schwarzschild and apply a cartesian coordinate transform
		name = 'black hole - Schwarzschild pseudocartesian',
		init = function(self, solver)
			solver.eqn:addGuiVars{
				{name = 'R', value = .002},	-- Schwarzschild radius
			}
		end,
		initState = function(self, solver)
			return [[
	const real R = gui_R;
	
	real r = real3_len(x);
	alpha = sqrt(1. - R/r);

	real3 xu = real3_scale(x, 1. / r);

	gamma_ll = sym3_add(
		sym3_ident(),
		sym3_scale(real3_outer(xu, xu), 1. / (r / R - 1)));
]]
		end,
	},
	{	-- Baumgarte & Shapiro, table 2.1, isotropic coordinates
		name = 'black hole - isotropic',
		init = function(self, solver)
			solver.eqn:addGuiVars{
				{name = 'R', value = .001},	-- Schwarzschild radius
			}
		end,
		initState = function(self, solver)
			solver:setBoundaryMethods'fixed'

			--momentum - isn't working
			local PU = {0,0,0}
			--local PU = {.1,0,0}
			--local PU = {1,0,0}
			
			--rotation
			--local JU = {0,0,0}
			--local JU = {0,0,.1}	//slowly comes to a stop and then forms a weird pattern and everything stops.
			local JU = {0,0,1}
			--local JU = {0,0,10}
	
			return template([[
	const real R = gui_R;
	real rSq = real3_lenSq(x);
	real r = sqrt(rSq);
	
	real alpha_num = 1. - .25 * R / r;
	real psi = 1. + .25 * R / r;		//sum of mass over radius
	real psi2 = psi * psi;
	real psi4 = psi2 * psi2;
	real psi8 = psi4 * psi4;

	alpha = alpha_num / psi;

	gamma_ll = sym3_scale(sym3_ident(), psi4);

	real3 PU = _real3(<?=clnumber(PU[1])?>, <?=clnumber(PU[2])?>, <?=clnumber(PU[3])?>);
	real3 JU = _real3(<?=clnumber(JU[1])?>, <?=clnumber(JU[2])?>, <?=clnumber(JU[3])?>);

	real rCubed = rSq * r;
	real r5 = rCubed * rSq;
	
	real3 lU = real3_scale(x, 1./r);

	//here's lowering it, using the diagonal metric specified above ...
	//gamma_ij = psi^4 delta_ij 
	//scaling by psi8 is the same as lowering two indexes

	//lower, only for cartesian:
	real3 lL = real3_scale(lU, psi4);
	real l_dot_P = real3_dot(lL, PU);

	sym3 ABar_boost_uu = sym3_scale(
		sym3_sub(
			sym3_add(real3_outer(PU, lU), real3_outer(lU, PU)),
			sym3_scale(sym3_sub(sym3_ident(), real3_outer(lU, lU)), l_dot_P))
		, 1.5/rSq);

	real3 rJU = real3_cross(JU, x);

	sym3 ABar_spin_uu = sym3_scale(
		sym3_add(real3_outer(x, rJU), real3_outer(rJU, x)),
		psi8 * 3. / r5);

	sym3 ABar_uu = sym3_add(ABar_boost_uu, ABar_spin_uu);

	//lower twice <-> scale by psi^8
	sym3 ABar_ll = sym3_scale(ABar_uu, psi8);

	//ATilde_ij = A_ij * exp(-4 phi)
	//K = 0, so A_ij = K_ij
	K_ll = sym3_scale(ABar_ll, 1. / psi2);

]], {
	JU = JU,
	PU = PU,
	clnumber = clnumber,
})
		end,
	},
	{
		name = 'binary black holes - isotropic',
		init = function(self, solver)
			solver.eqn:addGuiVars{
				{name = 'R1', value = .01},
				{name = 'R2', value = .01},
			}
		end,
		initState = function(self, solver)
			return [[
	const real R1 = gui_R1;
	const real R2 = gui_R2;
	
	real3 pos1 = x; pos1.x -= .25;
	real3 pos2 = x; pos2.x += .25;
	
	real r1 = real3_lenSq(pos1);
	real r2 = real3_lenSq(pos2);

	real psi = 1. + .25 * R1 / r1 + .25 * R2 / r2;
	real psi2 = psi * psi;
	real psi4 = psi2 * psi2;

	alpha = (1. - .25 * R1 / r1 - .25 * R2 / r2) / psi;

	gamma_ll = sym3_scale(sym3_ident(), psi4);
]]
		end,
	},
	{
		name = 'stellar model',
		init = function(self, solver)
			solver.eqn:addGuiVars{
				{name = 'bodyMass', value = .001},
				{name = 'bodyRadius', value = .1},
			}
		end,
		initState = function(self, solver)
			local bodies = table{
				{
					pos = {0,0,0},
					mass = solver.eqn.guiVars.bodyMass.value,
					radius = solver.eqn.guiVars.bodyRadius.value,
				}
			}
			
			return template([[
<? for _,body in ipairs(bodies) do
?>
	real3 pos = _real3(<?=table.map(body.pos, clnumber):concat', '?>);
	real3 ofs = real3_sub(x, pos);
	
	real r = real3_len(ofs);
	real3 l = real3_scale(ofs, 1./r);
	
	real m = r / gui_bodyRadius; m *= m * m; m = min(m, 1.); m *= gui_bodyMass;
	real R = 2. * m;

	alpha -= 2*m/r;
	gamma_ll = sym3_add(gamma_ll, sym3_scale(real3_outer(l,l), 1./(r/R - 1.)));

	if (r < gui_bodyRadius) {
		rho += gui_bodyMass / (4./3. * M_PI * gui_bodyRadius * gui_bodyRadius * gui_bodyRadius);
#if 0
		local r = math.sqrt(rSq)
		local rho0 = body.pressure or (body.mass / (4/3 * math.pi * body.radius * body.radius * body.radius))
		-- TOV pressure solution: http://physics.stackexchange.com/questions/69953/solving-the-tolman-oppenheimer-volkoff-tov-equation
		-- ... solution to constant pressure?  doesn't density decrease with radius? time to find a better source.
		pressure = pressure + (body.pressure or (rho0 * (
			(
				math.sqrt(1 - 2 * M / R) - math.sqrt(1 - 2 * M * r * r / (R * R * R))
			) / (
				math.sqrt(1 - 2 * M * r * r / (R * R * R)) - 3 * math.sqrt(1 - 2 * M / R)
			)
		)))
#endif
	}

<? end ?>

	alpha = sqrt(alpha);
]], {
	table = table,
	clnumber = clnumber,
	bodies = bodies,
})
		end,
	},
	--[=[
	{
		name = 'stellar model 3',
		getCodePrefix = function(self, solver)
			--[[
			earth radius = 6.37101e+6 m
			domain: 10x radius = 6.37101e+7 m
			earth mass = 5.9736e+24 kg = 5.9736e+24 * 6.6738480e-11 / 299792458^2 m
			earth mass = Em * G / c^2 in meters
			--]]
			local G = 6.6738480e-11 -- kg m^3/s^2
			local c = 299792458     -- m/s
			-- massInRadii is the order of 1e-9.  much more subtle than the default 1e-3 demo
			local earth = {radiusInM = 6.37101e+6, massInKg = 5.9736e+24}
			-- massInRadii is on the order of 1e-6
			local sun = {radiusInM = 6.960e+8, massInKg = 1.9891e+30}

			local planet = sun
			planet.massInM = planet.massInKg * G / c^2
			planet.massInRadii = planet.massInM / planet.radiusInM
			planet.radiusInCoords = .1
			planet.massInCoords = planet.massInRadii * planet.radiusInCoords
			for k,v in pairs(planet) do print(k,v) end

			local gridUnitsInM = planet.radiusInM / planet.radiusInCoords
			initNumRel{
				bodies={
					{pos = {0,0,0}, radius = planet.radiusInCoords, mass = planet.massInCoords},
	--]=]
	--[[
	2007 Alic et al "Efficient Implementation of finite volume methods in Numerical Relativity"
	1D Black Hole in wormhole form
	the paper says for >1D use the eqn 24 and do the 'stuffed black hole' fix
	the paper describes the metric in terms of eta, and then gives a relationship between r and eta ...
	but never labels their graphs, whether they are in terms of eta or r ...
	I tried eta and got an initial value of a lapse that looked, well, like a tanh would look (that is the lapse, according to the metric)
	 which is not what their graphs look like.  their graphs look like a sigmoid/tanh that goes from .9 to 1 ... not the range of tanh(eta) which goes from 0 to 1 on a domain of 0 to 8
	--]]
	{
		name = '1D black hole - wormhole form',
		init = function(self, solver)
			solver.eqn:addGuiVar{name='m', value=1}
		end,
		getCodePrefix = function(self, solver, getCodes)
			local m = solver.eqn.guiVars.m.value
		
			-- TODO refreshCodePrefix ... except this is called by that ...
			-- hmm, account for initState's defining the bounds
			solver.mins = require 'vec.vec3'(0,0,0) 
			solver.maxs = require 'vec.vec3'(8,8,8) 
			
			-- coords: t, eta, Omega
			-- g_tt = -(tanh eta)^2 = -alpha^2 <=> alpha = tanh eta
			-- g_eta_eta = 4 m^2 cosh(eta/2)^4
			-- g_Omega_Omega = 4 m^2 cosh(eta/2)^4
			-- using r = m/2 exp(Eta)
			-- m = mass
			-- TODO how should the finite volume scheme be modified to incorporate the volume element -- especially a dynamic volume element
			--local xNames = table{'eta', 'theta', 'phi'}
			-- TODO mapping to isotropic coordinates?
			-- or separate models or automatic transforms between isotropic and spherical?
			local xs = xNames:map(function(x) return symmath.var(x) end)
			
			-- [[ using eta as a coord
			local eta, theta, phi = xs:unpack()
			local alpha = symmath.tanh(eta)
			local g_eta_eta = 4 * m^2 * symmath.cosh(eta/2)^4
			local gamma = {g_eta_eta, 0, 0, g_eta_eta, 0, g_eta_eta}  
			--]]
			--[[
			local x,y,z = xs:unpack()
			local r = (x^2+y^2+z^2)^.5
			-- r = m/2 exp(eta)
			-- eta = log(2r/m)
			-- alpha = tanh(eta) = tanh(log(2r/m)) = sinh(eta) / cosh(eta) = (exp(eta) - exp(-eta)) / (exp(eta) + exp(-eta))
			-- = (2r/m - m/2r) / (2r/m + m/2r)
			-- = (4r^2 - m^2) / (4r^2 + m^2)
			local tanh_eta = (4 * r^2 - m^2) / (4 * r^2 + m^2)
			local alpha = tanh_eta
			-- or maybe the paper just initializes a lapse of 1? 
			--local alpha = 1
			-- 4 m^2 cosh(eta/2)^4 deta^2 = 4 m^2 cosh(log(2r/m)/2)^4 /r^2 dr^2
			-- = 4 (m/r)^2 (sqrt(2r/m) + sqrt(m/2r) )^4 dr^2
			-- but now we have an infinity as r->0 ...
			local g_rr = 4 * (m/r)^2 * (symmath.sqrt(2*r/m) + symmath.sqrt(.5*m/r))^4
			-- ... unless that "cosh eta / 2" means "cosh(eta)/2" .... but then why not just write 4/16 = 1/4 out front?  grrr this paper  ....
			-- 1/4 (m/r)^2 (2r/m + m/2r)^4
			--local g_rr = (m/r)^2 * (2*r/m + m/(2*r))^4 / 4
			local gamma = {g_rr, 0, 0, g_rr, 0, g_rr}
			--]]
			
			return initNumRel{
				solver = solver,
				getCodes = getCodes,
				vars = xs,
				alpha = alpha,
				gamma = gamma,
				K = {0,0,0,0,0,0},
			}
		end,
	},
--[[
2003 Alcubierre et al "Toward standard testbeds for numerical relativity"
2004 Bona et al "A symmetry-breaking mechanism..."
2009 Bona, Bona-Casas "Gowdy waves..." https://arxiv.org/pdf/0911.1208v1.pdf
Gaudy wave tests
1/sqrt(t) exp(Q/2) (-dt^2 + dz^2) + t (exp(P) dx^2 + exp(-P) dy^2)
Q & P are functions of t & z, only periodic in z
t = t0 exp(-tau/tau0)
dt = -t0/tau0 exp(-tau/tau0) dtau

specifically for this problem:
P = J0(2 pi t) cos(2 pi z)
Q = pi J0(2 pi) J1(2 pi) - 2 pi t J0(2 pi t) J1(2 pi t) cos(2 pi z)^2
		+ 2 pi^2 t^2 (J0(2 pi t)^2 + J1(2 pi t)^2 - J0(2 pi)^2 - J1(2 pi)^2)
2 pi t0 is the 20th root of the Bessel function J0 <=> t0 ~ 9.88

at time tau=0 <=> t=t0, and picking 2 pi t0 to be a root of J0 ...
Q = pi J0(2 pi) J1(2 pi) - 2 pi^2 t0^2 (J0(2 pi)^2 + J1(2 pi)^2)
--]]
	{
		name = 'Gowdy waves',
		getCodePrefix = function(self, solver)
			--solver.mins = require 'vec.vec3'(0,0,0)
			--solver.maxs = require 'vec.vec3'(10,10,10)
		end,
	},

-- 2003 Alcubierre et al "Toward standard testbeds for numerical relativity"
-- grids are from [-.5, .5] with 50 rho cells, rho = 1,2,4, and timesteps of .01/rho
	{
		name = 'testbed - robust',
		-- pick epsilon so epsilon^2 = 0
		-- eqn 4.1: pick epsilon from -1e-10 / rho^2 to 1e=10 / rho^2
		init = function(self, solver)
			solver.eqn:addGuiVar{name='epsilon', value=1e-10}
		end,
		resetState = function(self, solver)
			local epsilon = solver.eqn.guiVars.epsilon.value
			solver.eqn:fillRandom(epsilon)
		end,
	},
	{
		name = 'testbed - gauge wave',
		init = function(self, solver)
--			solver.mins = vec3(-.5, -.5, -.5)
--			solver.maxs = vec3(-.5, -.5, -.5)
			solver:setBoundaryMethods'periodic'
			solver.eqn:addGuiVars{
				{name='A', value=.1},	-- .1, .01
				{name='d', value=1},
			}
--			solver.eqn.guiVars.f.value = solver.eqn.guiVars.f.options:find'1'	-- set f=1
		end,
		initState = function(self, solver)
			return [[
	const real t = 0.;
	real theta = 2. * M_PI / gui_d * (x.x - t);
	real H = 1. + gui_A * sin(theta);
	alpha = sqrt(H);
	gamma_ll.xx = H;
	K_ll.xx = -M_PI * gui_A / gui_d * cos(theta) / alpha;
]]
		end,
	},
	{
		name = 'testbed - gauge wave - diagonal',
		init = function(self, solver)
			error"finishme"
		end,
	},
	{
		name = 'testbed - linear wave',
		init = function(self, solver)
-- TODO changing the range upon init causes something to freeze up ...
--			solver.mins = vec3(-.5, -.5, -.5)
--			solver.maxs = vec3(-.5, -.5, -.5)
			solver:setBoundaryMethods'periodic'
			solver.eqn:addGuiVars{
				{name='A', value=1e-8},
				{name='d', value=1},
			}
		end,
		initState = function(self, solver)
			return [[
	const real t = 0.;
	real theta = 2. * M_PI / gui_d * (x.x - t);
	real b = gui_A * sin(theta);
	gamma_ll.yy += b;
	gamma_ll.zz -= b;
	real db_dt = -2. * M_PI * gui_A / gui_d * cos(theta);
	K_ll.yy = .5 * db_dt;
	K_ll.zz = -.5 * db_dt;
]]
		end,
	},
	{
		name = 'testbed - linear wave - diagonal',
		init = function(self, solver)
			error"finishme"
		end,
	},
	{
		name = 'testbed - Gowdy',
		init = function(self, solver)
		
		end,
		getCodePrefix = function(self, solver, getCodes)
			local frac = symmath.frac
			local xs = xNames:map(function(x) return symmath.var(x) end)
			local x,y,z = xs:unpack()
			local t = 0
			--[[ general solution
			local lambda = -2 * pi * t * J0(2 * pi * t) * J1(2 * pi * t) * symmath.cos(2 * pi * z)^2
				+ 2 * pi^2 * t^2 * (J0(2 * pi * t)^2 + J1(2 * pi * t)^2)
				- 1/2 * ( (2 * pi)^2 * (J0(2 * pi)^2 + J1(2 * pi)^2) - 2 * pi * J0(2 * pi) * J1(2 * pi))
			--]]
			local lambda = -.5 * ( (2 * math.pi)^2 * (J0(2 * math.pi)^2 + J1(2 * math.pi)^2) - 2 * math.pi * J0(2 * math.pi) * J1(2 * math.pi))
			local alpha = t^symmath.frac(-1,4) * symmath.exp(lambda/4)
			local P = J0(2 * pi * t) * symmath.cos(2 * pi * z)
			local gamma = {t * symmath.exp(P), 0, 0, t * symmath.exp(-P), 0, 1/symmath.sqrt(t) * symmath.exp(lambda/2)}
			local K = {-1/2 * t^(1/4) * symmath.exp(P-lambda/4) * (1 + t * dP_dt), 0, 0, -1/2 * t^(1/4) * symmath.exp(-P-lambda/4) * (1 - t * dP_dt), 0, 1/4*t^(-1/4) * symmath.exp(lambda/4) * (1/t - dlambda_dt)}
			return initNumRel{
				solver = solver,
				getCodes = getCodes,
				vars = xs,
				alpha = alpha,
				gamma = gamma,
				K = K,
			}
		end,
	},
}:map(function(cl)
	return class(NumRelInitCond, cl)
end)

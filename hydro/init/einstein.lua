local class = require 'ext.class'
local table = require 'ext.table'
local symmath = require 'symmath'
local template = require 'template'
local vec3d = require 'vec-ffi.vec3d'
local clnumber = require 'cl.obj.number'
local InitCond = require 'hydro.init.init'
local constants = require 'hydro.constants'

local common = require 'hydro.common'
local xNames = common.xNames

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
		local code = symmath.export.C:toCode{
			output = vars,
		}, name
	
		-- ugh...
		code = code:gsub('sqrt%(', '(real)sqrt((real)')
		
		print(code)
		return code
	end
	return table.map(expr, function(v,k) return compileC(v,k,vars) end), name
end

local EinsteinInitCond = class(InitCond)

function EinsteinInitCond:init(args, ...)
	-- initialize analytical metric components to flat spacetime
	-- you still have to set 'initAnalytical=true'
	local Tensor = require 'symmath'.Tensor
	self.alpha0 = 1
	self.beta0_u = Tensor'^i'
	self.gamma0_ll = Tensor.metric().metric
	self.K0_ll = Tensor'_ij'
	
	--EinsteinInitCond.super.init(self, args, ...)
end

-- this is used atm
function EinsteinInitCond:getCodePrefix(solver)
	-- looks like all EFE solvers might need this
	-- maybe I should put it in InitCond?
	
	local alphaVar = symmath.var'alpha'
	-- TODO each eqn must be stated as a guiVar of the solver
	-- make a parent class or something for all Einstein field eqns
	local fGuiVar = solver.eqn.guiVars.f_eqn
	local fLuaCode = fGuiVar.options[fGuiVar.value]
	
	local f = assert(loadstring([[
local alpha, symmath = ...
local log = symmath.log
return ]]..fLuaCode))(alphaVar, symmath)
	f = symmath.clone(f)
	local dalpha_f = f:diff(alphaVar)()

	local codes = table()
	codes.f = compileC(f, 'f', {alphaVar})
	codes.f_alpha = compileC((f * alphaVar)(), 'f_alpha', {alphaVar})
	codes.f_alphaSq = compileC((f * alphaVar^2)(), 'f_alphaSq', {alphaVar})
	codes.dalpha_f = compileC(dalpha_f, 'dalpha_f', {alphaVar})
	codes.alpha_dalpha_f = compileC((alphaVar * dalpha_f)(), 'alpha_dalpha_f', {alphaVar})
	codes.alphaSq_dalpha_f = compileC((alphaVar^2 * dalpha_f)(), 'alphaSq_dalpha_f', {alphaVar})

	return codes:map(function(code,name,t)
		return 'real calc_'..name..'(real alpha) {\n\t'..code..'\n\treturn out1;\n}', #t+1
	end):concat'\n'		
end

-- TODO this is not in use yet
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
local function initEinstein(args)	
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
	-- make a parent class or something for all Einstein field eqns
	local fGuiVar = self.guiVars.f
	local fLuaCode = fGuiVar.options[fGuiVar.value]
	
	local f = assert(loadstring([[
local alpha, symmath = ...
local log = symmath.log
return ]]..fLuaCode))(alphaVar, symmath)
	f = symmath.clone(f)
	local dalpha_f = f:diff(alphaVar)()

	codes.f = compileC(f, 'f', {alphaVar})
	codes.dalpha_f = compileC(dalpha_f, 'dalpha_f', {alphaVar})

error'here'
	return table.map(codes, function(code,name,t)
		return 'real calc_'..name..'(real alpha) {\n\t'..code..'\n\treturn out1;\n}', #t+1
	end):concat'\n'
end

function EinsteinInitCond:buildFCCode(solver, diff)
	local alphaVar = symmath.var'alpha'
	local fGuiVar = solver.eqn.guiVars.f_eqn
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
	{
		name = 'Minkowski',
		-- flag for determining whether to initialize variables (esp the derivative variables) from analytical expressions, or whether to use finite difference via initDerivs
		-- don't use initAnalytical, since not all eqn.einstein subclasses use it
		--initAnalytical = true,
		getInitCondCode = function(self, solver)
			return [[
	alpha = 1.;
	beta_u = real3_zero;
	gamma_ll = coord_g_ll(x);
	K_ll = sym3_zero;
]]
		end,
	},

	-- 2010 Baumgarte, Shapiro "Numerical Relativity ...", section 9.1.2
	-- 1997 Alcubierre "The appearance of coorindate shocks in hyperbolic formalisms of General Relativity".
	{
		name = 'gaussian perturbation',
		--[[
		args:
			center (in grid coordinates)
			H
			sigma
		--]]
		init = function(self, solver, args)
			self.initCondArgs = args
			EinsteinInitCond.init(self, solver, args)
		end,
		createInitStruct = function(self, solver)
			EinsteinInitCond.createInitStruct(self, solver)
			local args = self.initCondArgs

			local size = solver.maxs.x - solver.mins.x
			self.center = args and args.center or {0,0,0}
			
			-- 1997 Alcubierre uses amplitude of 5 on a grid size of 300 with dx=1
			--local H = 5 / 300 * size
			--local sigma = 10 / 300 * size
			-- 1998 Bona et al use amplitude of .05 on a grid size of 100 with dx=.01
			--local H = .05 * size
			--local sigma = math.sqrt(.05) * size
			-- messing with it on my own
			local H = args and args.H or .01 * size
			local sigma = args and args.sigma or .1 * size

			self:addGuiVars{
				{name = 'H', value = H},
				{name = 'sigma', value = sigma},
			}
		end,
		getInitCondCode = function(self, solver)
			-- this has to make use of the coordinate metric
			-- solver.coord.g
			
			return template([[
	real3 center = coordMap(_real3(<?=clnumber(initCond.center[1])
								?>, <?=clnumber(initCond.center[2])
								?>, <?=clnumber(initCond.center[3])?>));
	real3 d = real3_sub(xc, center);
	real s = real3_lenSq(d);

	const real H = initCond->H;
	const real sigma = initCond->sigma;
	const real sigma2 = sigma * sigma;
	const real sigma4 = sigma2 * sigma2;
	real h = H * exp(-s / sigma2);

	//h,i = (H exp(-(x-xc)^2 / sigma^2)),i
	// = -2 h / sigma^2 (x_i-c_i)
	real3 dh = real3_real_mul(real3_sub(x, mids), -2. * h / sigma2);

	sym3 delta_ll = sym3_ident;

	//h,ij = h,i,j = (-2 h / sigma^2 (x_i-c_i)),j
	// = -2 h,j / sigma^2 (x_i-c_i) - 2 h delta_ij / sigma^2
	// = -2 (-2 h / sigma^2 (x-c)_i (x-c)_j) / sigma^2 - 2 h delta_ij / sigma^2
	// = 4 h (x-c)_i (x-c)_j / sigma^4 - 2 h delta_ij / sigma^2
	sym3 d2h = sym3_sub(
		sym3_real_mul(real3_outer(d), 4. * h / sigma4),
		sym3_real_mul(sym3_ident, 2. * h / sigma2));

#if 0
	alpha = .5 * (
		(1 + sqrt(1 + kappa)) * sqrt((1-dh.x)/(1+dh.x))
		- kappa / (1 + sqrt(1 + kappa)) * sqrt((1+dh.x)/(1-dh.x))
	);
#endif

	gamma_ll = sym3_sub(delta_ll, real3_outer(dh));
	K_ll = sym3_real_mul(d2h, -1./sqrt(1. - real3_lenSq(dh)));

//plane wave vs bubble perturbation
//TODO make these two separate init conds
//enable this if you want the ADM 3D run in 1D to match the ADM 1D's
//disable this if you want things to run in higher dimensions
#if	dim == 1
	gamma_ll = _sym3(gamma_ll.xx, 0,0,1,0,1);
	K_ll = _sym3(K_ll.xx, 0,0,0,0,0);
#endif
]],			{
				initCond = self,
				clnumber = require 'cl.obj.number',
			})
		end,
	},
	-- from 2012 Alic et. al. "Conformal and covariant formulations of the Z4 system with constraint-violation damping"
	{
		name = 'plane gauge wave',
		guiVars = {
			{name = 'A', value = .1},
			{name = 'L', value = 1},
		},
		getInitCondCode = function(self, solver)
			return [[
	real h = 1. - initCond->A * sin((2. * M_PI / initCond->L) * x.x);
	alpha = sqrt(h);
	gamma_ll.xx = h;

]]
		end,
	},
	
	-- from 2011 Alcubierre, Mendez "Formulations of the 3+1 evolution equations in curvilinear coordinates"
	-- Appendix A, eqns 7.1 - 7.5
	{
		name = 'pure gauge wave',
		initAnalytical = true,
		guiVars = {
			{name = 'alpha0', value = .01},
			{name = 'r0', value = 5},
			{name = 'sigma', value = 1},
		},
		init = function(self, solver, args)
			EinsteinInitCond.init(self, solver, args)
			
			local coord = solver.coord
			assert(coord)

			local symmath = require 'symmath'
			local var = symmath.var
			local exp = symmath.exp
			
			-- TODO make a separate OpenCL initCond_t structure
			-- this is the alpha0 of the init state params -- wave amplitude -- not the init cond alpha0 used for the initial alpha value
			local param_alpha0 = var'initCond->alpha0'
			local param_r0 = var'initCond->r0'
			local param_sigma = var'initCond->sigma'
	
			-- variables used by initAnalytical
			local r = coord.vars.r 
			local rplus = (r + param_r0) / param_sigma
			local rminus = (r - param_r0) / param_sigma
			local gminus = exp(-rminus^2)
			local gplus = exp(-rplus^2)
			
			-- radial (/ first dimension) perturbation of alpha
			self.alpha0 = 1 + param_alpha0 * r^2 / (1 + r^2) * (gplus + gminus)
		end,
	},

	{
		name = 'scalar field',
		guiVars = {
			{name = 'alpha0', value = .01},
			{name = 'r0', value = 5},
			{name = 'sigma', value = 1},
		},
		getInitCondCode = function(self, solver)
			assert(solver.eqn.useScalarField, "you need to enable 'useScalarField' in the eqn ctor args")
			
			local symmath = require 'symmath'
			local var = symmath.var
			local Tensor = symmath.Tensor

			local r = solver.coord.vars.r
			local r0 = var'initCond->r0'
			local sigma = var'initCond->sigma'
			local alpha0 = var'initCond->alpha0'
			local rplus = (r + r0) / sigma
			local rminus = (r - r0) / sigma
			local gminus = symmath.exp(-rminus * rminus)
			local gplus = symmath.exp(-rplus * rplus)
			local Phi = alpha0 * r^2 / (1 + r^2) * (gplus + gminus)
			
			local baseCoords = solver.coord.baseCoords
			local Psi = Tensor('_i', function(i)
				return Phi:diff(baseCoords[i])()
			end)
			
			-- Pi is 1/alpha (Phi_,t - beta^i Psi_i)
			local Phi_t = 0
			local alpha = var'alpha'	-- ADM var
			local beta = Tensor('^i', function(i) 
				return var('beta_u.s'..(i-1))
			end)
			local Pi = 1/alpha * (Phi_t - beta'^i' * Psi'_i')()

			local compile = function(...)
				local code = solver.eqn:compile(...)
				code = code:gsub('pt%.', 'x.')
				return code
			end

		return template([[
	Phi = cplx_from_real(<?=compile(Phi)?>);
	real3 re_Psi_l = _real3(
		<?=compile(Psi[1])?>,
		<?=compile(Psi[2])?>,
		<?=compile(Psi[3])?>);
	Psi_l = cplx3_from_real3(re_Psi_l);
	Phi = cplx_from_real(<?=compile(Phi)?>);
]], 	{
			Phi = Phi,
			Psi = Psi,
			Pi = Pi,
			symmath = symmath,
			compile = compile,
		})
		end,
	},


	{
		name = 'Alcubierre warp bubble',
		--[[
		args:
			R = warp bubble radius
			sigma = warp bubble thickness
			speed = warp bubble speed
		--]]
		init = function(self, solver, args)
			self.initCondArgs = args
			EinsteinInitCond.init(self, solver, args)
		end,
		createInitStruct = function(self, solver)
			EinsteinInitCond.createInitStruct(self, solver)
			local args = self.initCondArgs

			self:addGuiVars{
				{name = 'R', value = args and args.R or .5},
				{name = 'sigma', value = args and args.sigma or 8},
				{name = 'speed', value = args and args.speed or .1},
			}
		end,
		getInitCondCode = function(self, solver)
			return [[
	real x_s = 0;	//speed * t
	real v_s = initCond->speed;
	
	real3 y = xc; y.x -= x_s;
	real r_s = real3_len(y);

#define cosh(x)		(.5 * (exp(x) + exp(-x)))
#define dtanh(x) 	(1./(cosh(x)*cosh(x)))

	real fnum = tanh(initCond->sigma * (r_s + initCond->R)) - tanh(initCond->sigma * (r_s - initCond->R));
	real fdenom = 2 * tanh(initCond->sigma * initCond->R);
	real f = fnum / fdenom;

	beta_u.x = -v_s * f;

	real3 dx_r_s = real3_real_mul(y, 1. / r_s);

	real3 dx_f = real3_real_mul(dx_r_s, 
		initCond->sigma * (
			dtanh(initCond->sigma * (r_s + initCond->R)) 
			- dtanh(initCond->sigma * (r_s - initCond->R))
		) / fdenom);

	alpha = 1;

	K_ll.xx = -v_s * dx_f.x / alpha;
	K_ll.xy = -v_s * dx_f.y / (2 * alpha);
	K_ll.xz = -v_s * dx_f.z / (2 * alpha);
]]
		end,
	},
	
	
	
	{	
		name = 'black hole - Schwarzschild',
		init = function(self, solver, args)
			self.initCondArgs = args
			EinsteinInitCond.init(self, solver, args)
		end,
		createInitStruct = function(self, solver)
			EinsteinInitCond.createInitStruct(self, solver)
			args = self.initCondArgs or {}
				
			-- TODO bodies
			self:addGuiVars{
				{name = 'R', value = args and args.R or 1},
				{name = 'x', value = args and args.x or 0},
				{name = 'y', value = args and args.y or 0},
				{name = 'z', value = args and args.z or 0},
			}
		end,
		getInitCondCode = function(self, solver)
			return template([[
	const real R = initCond->R;
	real3 center = _real3(
		initCond->x,
		initCond->y,
		initCond->z);
	real3 xrel = real3_sub(xc, center);

	real r = real3_len(xrel);

	//real one_minus_R_over_r = (fabs(r) < 1e-3 || fabs(r - R) < 1e-3) ? 1. : 1. - R / r;
	real one_minus_R_over_r = 1. - R / r;
	alpha = sqrt(one_minus_R_over_r);

	//bssnok-fd init cond expects the metric in non-coord normalized basis
	//...and all the finite-volume solvers (adm3d, etc) only work with cartesian grids
	gamma_ll = sym3_zero;
	gamma_ll.xx = 1. / one_minus_R_over_r;
	gamma_ll.yy = 1. / one_minus_R_over_r;
	gamma_ll.zz = 1. / one_minus_R_over_r;
]], 	{
			solver = solver,
		})
		end,
	},
	
	
	
	
	{	-- Baumgarte & Shapiro, table 2.1, isotropic coordinates
		-- also looking at ch 12.2 on binary black hole initial data
		--[[
		this is going to require one extra computation, similar to initDerivs,
		for performing the inverse laplacian to compute the 'u' parameter
		which is used to compute the psi parameter
		--]]
		name = 'black hole - isotropic',
		useBSSNVars = true,	-- getInitCondCode writes to BSSN vars
		--[[
		args:
			bodies:
				R = Schwarzschild radius
				P_u = linear momentum 
				S_u = angular momentum 
				pos = separation
		--]]
		init = function(self, solver, args)
			EinsteinInitCond.init(self, solver, args)
			
			args = args or {}
			
			local v = self.guiVars
			self.bodies = args.bodies or {
				{
					R = .0001,
					P_u = {0,0,0},
					S_u = {0,0,.01},
					pos = {0,0,0},
				},
			}
			
			-- TODO allow rescaling grid within initial conditions.
			--[[ TODO use these instead of bodies ... and have a parameter for # of bodies ... and add/subtract bodies ... 
			for i,body in ipairs(self.bodies) do
				local bodyPrefix = i..'.'
				solver.eqn:addGuiVar{name = bodyPrefix..'R', value = body.R}
				for j,xj in ipairs(xNames) do
					solver.eqn:addGuiVar{name = bodyPrefix..'P'..xj, value = body.P_u[j]}
				end
				for j,xj in ipairs(xNames) do
					solver.eqn:addGuiVar{name = bodyPrefix..'S'..xj, value = body.S_u[j]}
				end
				for j,xj in ipairs(xNames) do
					solver.eqn:addGuiVar{name = bodyPrefix..'pos'..xj, value = body.S_u[j]}
				end
				solver.eqn:addGuiVar{name = bodyPrefix..'dist', value = body.dist}
			end
			--]]
		end,
		getInitCondCode = function(self, solver)
			--solver:setBoundaryMethods'fixed'
			--solver:setBoundaryMethods'linear'

			--[[
			I'm trying to follow 1997 Brandt & Brugmann, but here's what I've gathered:
			1) calculate K_PS^ij for each puncture.  two lowerings are involved in computing the spatial cross product.  these require phi -- which hasn't been calculated yet.  (should I pretend they are lowered with an identity metric?)
			2) use the sum of K_PS^ij's to compute K^ij at each point
			3) compute alpha = 1 / (.5 sum m_i / (r - r_i))
			4) compute beta = 1/8 alpha^7 K^ab K_ab ... which again requires phi, which is what we're trying to solve for ... (is this a gradient descent or a predictor-corrector algorithm?) 
			5) compute (converge?) u = lap^-1 beta / (1 + alpha u)^7
			6) phi = 1/alpha + u
			7) and then there's something about the spherical coordinate inversion ...
			
			Here's the answer to the circular definitions of index raising/lowering with a metric that uses the variables we are solving for ...
			B&S 12.2.1: "...so that gammaTilde_ij = eta_ij..."
			So all raising and lowering throughout these computations uses eta_ij?
			This statement makes sense to consider raising/lowering of ABar and gammaTilde using eta_ij.
			However it does not explain the raising/lowering of non-bar terms like n^i, S^i, P^i etc
			And I would simply change the ABar computations from lower to upper, except even they include a term n_i S^i which will still require a metric.

			... or (for beta) just use 2015 Baumgarte's advice and initialize to beta=0 ... which I'm doing anyways

			--]]
			return template([[
	//real sln_oneOverAlpha = 0.;
	
	real psi = 1.;
	real alpha_num = 1.;
	<? for _,body in ipairs(bodies) do ?>{
		const real R = <?=clnumber(body.R)?>;
		real3 pos = _real3(<?=table(body.pos):mapi(clnumber):concat', '?>);
		real3 xrel = real3_sub(xc, pos);
		real r = real3_len(xrel);
		
		// the correct binary black hole way: B&S 12.51
		//sln_oneOverAlpha += .5 * <?=clnumber(body.R)?> / r;
		
		//the time-symmetric way:
		psi += .25 * R / r;
		alpha_num -= .25 * R / r;
	
	}<? end ?>		//sum of mass over radius

	//correct binary black hole way
	//psi = 1 + 1/alpha + u	-- B&S 12.50

	//next solve u via 
	// Baumgarte & Shapiro add an extra 1 to psi's calculations as compared to the original Brandt & Brugmann equation
	//DBar^2 u = -1/8 alpha^7 ABar_ij ABar^ij / ( alpha (1 + u) + 1)^7 	<-- Baumgarte & Shapiro eqns 12.52, 12.53
	//...factoring out the alpha ...
	//DBar^2 u = -1/8 ABar_ij ABar^ij / (1 + 1/alpha + u)^7

	//time-symmmetric way:
	alpha = alpha_num / psi;
	real psi2 = psi * psi;
	real psi4 = psi2 * psi2;
	
	W = psi2;

	<? for _,body in ipairs(bodies) do ?>{

		real3 pos = _real3(<?=table(body.pos):mapi(clnumber):concat', '?>);
		real3 xrel = real3_sub(xc, pos);
		real r = real3_len(xrel);
		real rSq = r * r;
		real rCubed = rSq * r;

		// upper is cartesian coordinates
		// metric is isotropic
		real3 n_u = r == 0 ? _real3(0,0,1) : real3_real_mul(xrel, 1./r);
		real3 n_l = real3_real_mul(n_u, psi4);

		real3 P_u = _real3(<?=table(body.P_u):mapi(clnumber):concat', '?>);
		real3 S_u = _real3(<?=table(body.S_u):mapi(clnumber):concat', '?>);

		real3 P_l = real3_real_mul(P_u, psi4);
		
		real n_dot_P = real3_dot(n_l, P_u);
		
		//Alcubierre 3.4.22
		//Bowen-York extrinsic curvature

		sym3 ABar_boost_ll = sym3_real_mul(
			sym3_sub(
				sym3_real_mul(sym3_from_real3x3(real3_real3_outer(P_l, n_l)), 2.),
				sym3_real_mul(sym3_sub(real3_outer(n_l), sym3_ident), n_dot_P)
			), 1.5 / rSq);

		//Levi-Civita density is det gamma for conformal metric, whose det is 1, so a cross product works with covariant Levi-Civita
		real3 S_cross_n_l = real3_cross(S_u, n_u);

		sym3 ABar_spin_ll = sym3_real_mul(
			sym3_real_mul(sym3_from_real3x3(real3_real3_outer(S_cross_n_l, n_l)), 2.),
			3. / rCubed);

		ABar_LL = sym3_add(
			ABar_LL, 
			sym3_rescaleFromCoord_LL(
				sym3_add(ABar_boost_ll, ABar_spin_ll),
				x)
		);
	}<? end ?>

]], 		{
				bodies = self.bodies,
				clnumber = clnumber,
			})
		end,
	},

	{
		-- How to provide a single initial condition for any coordinate system?
		-- I could always do coordinate transforms based on e_i^I, however this wouldn't take into account change-of-coordinates that are required for things like spherical <-> isotropic pseudo-cartesian
		-- In the mean time I'll make a new init cond for each geometry type.
		-- Spherical and sphere-log-polar should still be interoperable.
		name = 'black hole - Schwarzschild - spherical',
		initAnalytical = true,
		guiVars = {
			{name = 'rs', value = .1},
		},
		init = function(self, solver, args)
			EinsteinInitCond.init(self, solver, args)
			
			local symmath = require 'symmath'
			setfenv(1, setmetatable({}, {
				__index = function(t,k)
					local v = symmath[k] if v ~= nil then return v end
					local v = _G[k] if v ~= nil then return v end
				end,
			}))
			
			local r = solver.coord.vars.r
			local theta = solver.coord.baseCoords[2]	-- only true for sphere and sphere-log-polar
			local rs = var'initCond->rs'
			
			self.alpha0 = sqrt(1 - rs/r)
			self.gamma0_ll = Tensor('_ij', 
				{1/(1 - rs/r), 0, 0},
				{0, r^2, 0},
				{0, 0, r^2 * sin(theta)^2}
			)
		end,
	},

	{
		name = 'black hole - Schwarzschild isotropic - spherical',
		initAnalytical = true,
		guiVars = {
			{name = 'rs', value = .1},
		},
		init = function(self, solver, args)
			EinsteinInitCond.init(self, solver, args)
			
			local symmath = require 'symmath'
			setfenv(1, setmetatable({}, {
				__index = function(t,k)
					local v = symmath[k] if v ~= nil then return v end
					local v = _G[k] if v ~= nil then return v end
				end,
			}))
			
			local r = solver.coord.vars.r
			local theta = solver.coord.baseCoords[2]	-- only true for sphere and sphere-log-polar
			local rs = var'initCond->rs'
			
			self.alpha0 = sqrt( (1 - rs/(4*r)) / (1 + rs/(4*r)) )
			self.gamma0_ll = (Tensor('_ij', 
				{1, 0, 0},
				{0, r^2, 0},
				{0, 0, r^2 * sin(theta)^2}
			) * (1 + rs/(4*r))^4)()
		end,
	},


	-- SENR's BoostedSchwarzschild initial data
	{
		name = 'black hole - boosted Schwarzschild',
		init = function(self, solver)
			EinsteinInitCond.init(self, solver, args)
		end,
		getInitCondCode = function(self, solver)
			return [[
	{
		const real vz = .1;
		const real M = 1.;
		const real xB = 0.;
		const real yB = 0.;
		const real zB = 0.;
		real vzsq = vz * vz;
		real vz_toThe4 = vzsq * vzsq;
		real LF = 1. / sqrt(1. - vzsq);
		real LFsq = LF*LF;
		real LF_toThe3 = LFsq * LF;
		real rB = sqrt(xB*xB + yB*yB + LFsq*zB*zB);
		real psiB = 1. + M / (2 * rB);
		real psiBsq = psiB * psiB;
		real psiB_toThe4 = psiBsq * psiBsq;
		real M_plus_2_rB = M + 2*rB;
		real M_plus_2_rB_sq = M_plus_2_rB * M_plus_2_rB;
		real M_plus_2_rB_toThe3 = M_plus_2_rB_sq * M_plus_2_rB;
		real M_plus_2_rB_toThe5 = M_plus_2_rB_toThe3 * M_plus_2_rB_sq;
		real M_plus_2_rB_toThe6 = M_plus_2_rB_toThe3 * M_plus_2_rB_toThe3;
		real M_plus_2_rB_toThe7 = M_plus_2_rB_toThe5 * M_plus_2_rB_sq;
		real M_plus_2_rB_toThe9 = M_plus_2_rB_toThe6 * M_plus_2_rB_toThe3;
		real M_plus_2_rB_toThe12 = M_plus_2_rB_toThe6 * M_plus_2_rB_toThe6;
		real M_minus_2_rB = M - 2*rB;
		real M_minus_2_rB_sq = M_minus_2_rB * M_minus_2_rB;
		real M_minus_2_rB_toThe4 = M_minus_2_rB_sq * M_minus_2_rB_sq;
		real M_minus_4_rB = M - 4.*rB;
		real rBsq = rB * rB;
		real rB_toThe4 = rBsq * rBsq;
		real rB_toThe8 = rB_toThe4 * rB_toThe4;
		real BBsq = M_plus_2_rB_toThe6 - 16. * M_minus_2_rB_sq * rB_toThe4 * vzsq;
		real BB = sqrt(BBsq);
		real BB_toThe3 = BB * BBsq;
		real CC = M_plus_2_rB_toThe6 - 8. * M_minus_2_rB_sq * rB_toThe4 * vzsq;
		beta_u.z = -1./(BB*BB) * (
			M * vz 
			* (M*M + 6*M*rB + 16*rBsq) 
			* (M*M*M + 6*M*M*rB + 8*M*rB*rB + 16*rB*rBsq)
		);
		B_u.z = -(
			64*M*M_minus_4_rB*M_minus_2_rB*rBsq*M_plus_2_rB_toThe5 * vzsq*zB
		) / (
			M_plus_2_rB_toThe12 - 32*M_minus_2_rB_sq * rB_toThe4 * M_plus_2_rB_toThe6 * vzsq + 256*M_minus_2_rB_toThe4 * rB_toThe8 * vz_toThe4
		);
		real gtzz = LFsq * BB*BB / M_plus_2_rB_toThe6;
		real Atxx = (64*LF*M*M_minus_4_rB*rBsq*vz*CC*zB) / (3*M_plus_2_rB_toThe3*BB_toThe3);
		real Atxz = -(32*LF*M*M_minus_4_rB*rBsq*vz*xB) / (M_plus_2_rB_toThe3*BB);
		real Atyz = -(32*LF*M*M_minus_4_rB*rBsq*vz*yB) / (M_plus_2_rB_toThe3*BB);
		real Atzz = -(128*LF_toThe3*M*M_minus_4_rB*rBsq*vz*CC*zB) / (3*M_plus_2_rB_toThe9*BB);

		real r0 = x.x;
		real theta0 = x.y;
		real phi0 = x.z;
		real r0sq = r0 * r0;
		
		real cos_theta0 = cos(theta0);
		real cos_theta0_sq = cos_theta0 * cos_theta0;
		real sin_theta0 = sin(theta0);
		real sin_theta0_sq = sin_theta0 * sin_theta0;

		sym3 gPhys00DD;
		gPhys00DD.xx = psiB_toThe4 * .5 * (1. + gtzz + (-1. + gtzz) * cos(2. * theta0));
		gPhys00DD.xy = psiB_toThe4 * (1. - gtzz) * r0 * cos_theta0 * sin(theta0);
		gPhys00DD.xz = 0.;
		gPhys00DD.yy = psiB_toThe4 * .5 * r0sq * (1. + gtzz - (-1. + gtzz) * cos(2. * theta0));
		gPhys00DD.yz = 0.;
		gPhys00DD.zz = psiB_toThe4 * r0sq * sin(theta0) * sin(theta0);
		gamma_ll = gPhys00DD;

		sym3 APhys00DD;
		APhys00DD.xx = psiB_toThe4 * (
			.5*(Atxx + Atzz + (-Atxx + Atzz)*cos(2*theta0) + 2*(Atxz*cos(phi0) + Atyz*sin(phi0))*sin(2*theta0))
		);
		APhys00DD.xy = psiB_toThe4 * (
			r0*cos(2*theta0)*(Atxz*cos(phi0) + Atyz*sin(phi0)) + (Atxx - Atzz)*r0*cos_theta0*sin(theta0)
		);
		APhys00DD.xz = psiB_toThe4 * (
			r0*cos_theta0*(Atyz*cos(phi0) - Atxz*sin(phi0))*sin(theta0) 
		);
		APhys00DD.yy = psiB_toThe4 * (
			r0sq*(Atxx*cos_theta0_sq + Atzz*sin_theta0_sq - (Atxz*cos(phi0) + Atyz*sin(phi0))*sin(2*theta0))	
		);
		APhys00DD.yz = psiB_toThe4 * (
			r0sq*(-Atyz*cos(phi0) + Atxz*sin(phi0))*sin_theta0_sq
		);
		APhys00DD.zz = psiB_toThe4 * (
			Atxx*r0sq*sin_theta0_sq	
		);

		real K = (32*LF * M * vz * (M_plus_2_rB_toThe7 - 32*M_minus_2_rB_sq * (M - rB) * rB_toThe4 * vzsq)*rBsq * zB)
			/ (M_plus_2_rB_toThe3*BB_toThe3);
	
		real det_gamma = sym3_det(gamma_ll);
		real det_gammaBar = calc_det_gammaBar(x); 
		real exp_neg4phi = cbrt(det_gammaBar / det_gamma);

		sym3 A_ll = sym3_real_mul(APhys00DD, 1. / exp_neg4phi);
		K_ll = sym3_add(A_ll, sym3_real_mul(gamma_ll, K));
	}
]]
		end,
	},

	{	-- initial data from SENR/NumPy jupyter notebooks:
		-- https://hub.mybinder.org/user/zachetienne-nrpytutorial-a0u5aw7y/notebooks/Tutorial-ADM_Initial_Data-Brill-Lindquist.ipynb
		name = 'black hole - Brill Lindquist',
		init = function(self, solver)
			EinsteinInitCond.init(self, solver, args)
			
			args = args or {}
			local v = self.guiVars
			self.bodies = args.bodies or {
				{
					R = .0001,
					P_u = {0,0,0},
					S_u = {0,0,0},
					pos = {0,0,0},
				},
			}
		end,
		getInitCondCode = function(self, solver)
			return template([[
<? local clnumber = require 'cl.obj.number' ?>	
	real psi = 1.;
	<? for _,body in ipairs(bodies) do ?>{
		const real R = <?=clnumber(body.R)?>;
		real3 pos = _real3(<?=clnumber(body.pos[1])?>, <?=clnumber(body.pos[2])?>, <?=clnumber(body.pos[3])?>);
		real3 xrel = real3_sub(xc, pos);
		real r = real3_len(xrel);
		psi += .25 * R / r;
	}<? end ?>
	real psi2 = psi*psi;
	gamma_ll = sym3_real_mul(sym3_ident, psi2*psi2);
	alpha = 1. / psi2;
]], 		{
				bodies = self.bodies,
			})
		end,
	},

	{	-- another based on SENR
		-- TODO make a way to provide BSSN variables rather than ADM variables
		-- this way I don't have to worry about scaling and unscaling of the coordinate basis vectors
		getInitCondCode = function(self, solver)
			return [[
	real bScale = 1.;
	real bScale2 = bScale * bScale;
	real r0 = real3_len(xc);
	real theta0 = atan2(xc.y, xc.x);
	real r02 = r0 * r0;
	real cos_theta0 = cos(theta0);
	real M1 = 1.;
	real M2 = 1.;
	real psi0 = 1 + .5 * (
		M1 / sqrt(bScale2 + r02 - 2*bScale*r0*cos_theta0)
		+ M2 / sqrt(bScale2 + r02 + 2*bScale*r0*cos_theta0)
	);
	U->alpha = 1.;	//psi0?
	//psi^2 = exp(-2 phi) = W
	U->W = psi0 * psi0;
	U->K = 0.;
	U->beta_U = real3_zero;
	U->B_U = real3_zero;
	U->LambdaBar_U = real3_zero;
	U->epsilon_LL = sym3_zero; 
	U->ABar_LL = sym3_zero;
]]
		end,
	},
	
	-- based on ch.23 of 1973 Misner, Thorne, Wheeler "Gravitation"
	-- TODO add support for multiple bodies
	{
		name = 'stellar model',
		guiVars = {
			{name = 'bodyMass', value = .001},
			{name = 'bodyRadius', value = .1},
		},
		getInitCondCode = function(self, solver)
			local bodies = table{
				{
					pos = {0,0,0},
					mass = self.guiVars.bodyMass.value,
					radius = self.guiVars.bodyRadius.value,
				}
			}
			
			return template([[
<? for _,body in ipairs(bodies) do
?>
	real3 pos = _real3(<?=table.map(body.pos, clnumber):concat', '?>);
	real3 ofs = real3_sub(xc, pos);
	
	real r = real3_len(ofs);
	real3 l = real3_real_mul(ofs, 1./r);
	
	real m = r / initCond->bodyRadius; m *= m * m; m = min(m, 1.); m *= initCond->bodyMass;
	real R = 2. * m;

	alpha -= 2*m/r;
	gamma_ll = sym3_add(gamma_ll, sym3_real_mul(real3_outer(l), 1./(r/R - 1.)));

	if (r < initCond->bodyRadius) {
		rho += initCond->bodyMass / (4./3. * M_PI * initCond->bodyRadius * initCond->bodyRadius * initCond->bodyRadius);
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
]], 		{
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
			local G = constants.gravitationalConstant_in_m_per_s  -- kg m^3/s^2
			local c = constants.speedOfLight_in_m_per_s -- m/s
			-- massInRadii is the order of 1e-9.  much more subtle than the default 1e-3 demo
			local earth = {radiusInM = constants.EarthRadius_in_m, massInKg = constants.EarthMass_in_kg}
			-- massInRadii is on the order of 1e-6
			local sun = {radiusInM = constants.SolarRadius_in_m, massInKg = constants.SolarMass_in_kg}

			local planet = sun
			planet.massInM = planet.massInKg * G / c^2
			planet.massInRadii = planet.massInM / planet.radiusInM
			planet.radiusInCoords = .1
			planet.massInCoords = planet.massInRadii * planet.radiusInCoords
			for k,v in pairs(planet) do print(k,v) end

			local gridUnitsInM = planet.radiusInM / planet.radiusInCoords
			initEinstein{
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
			EinsteinInitCond.init(self, solver, args)
			
			solver.eqn:addGuiVar{name='m', value=1}
		end,
		getCodePrefix = function(self, solver)
			local m = self.guiVars.m.value
		
			-- TODO refreshCodePrefix ... except this is called by that ...
			-- hmm, account for initCond's defining the bounds
			solver.mins = vec3d(0,0,0) 
			solver.maxs = vec3d(8,8,8) 
			
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
			
			return initEinstein{
				solver = solver,
				-- looks like I am no longer passing this into this function...
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

TODO I now have a Bessel function routine in hydro/math.cl
--]]
	{
		name = 'Gowdy waves',
		getCodePrefix = function(self, solver)
			--solver.mins = vec3d(0,0,0)
			--solver.maxs = vec3d(10,10,10)
		end,
	},

-- 2003 Alcubierre et al "Toward standard testbeds for numerical relativity"
-- grids are from [-.5, .5] with 50 rho cells, rho = 1,2,4, and timesteps of .01/rho
	{
		name = 'testbed - robust',
		-- pick epsilon so epsilon^2 = 0
		-- eqn 4.1: pick epsilon from -1e-10 / rho^2 to 1e=10 / rho^2
		-- the SENR/NumPy uses an epsilon of .02
		init = function(self, solver)
			EinsteinInitCond.init(self, solver, args)
			
			solver.eqn:addGuiVar{name='epsilon', value=1e-10}
		end,
		getInitCondCode = function(self, solver)
			return [[
	alpha = 1. + initCond->epsilon * U.ptr[0];
	for (int j = 0; j < 6; ++j) {
		gamma_ll.s[j] = initCond->epsilon * U.ptr[j+1];
	}
	gamma_ll = sym3_add(gamma_ll, conn_g_ll(x));
	for (int j = 0; j < 6; ++j) {
		K_ll.s[j] = initCond->epsilon * U.ptr[j+7];
	}
	for (int j = 0; j < 3; ++j) {
		beta_u.s[j] = initCond->epsilon * U.ptr[j+10];
	}
]]
		end,
	},
	{
		name = 'testbed - gauge wave',
		createInitStruct = function(self, solver)
			EinsteinInitCond.createInitStruct(self, solver)

--			solver.mins = vec3d(-.5, -.5, -.5)
--			solver.maxs = vec3d(-.5, -.5, -.5)
			solver:setBoundaryMethods'periodic'
			self:addGuiVars{
				{name='A', value=.1},	-- .1, .01
				{name='d', value=1},
			}
--			self.guiVars.f.value = self.guiVars.f.options:find'1'	-- set f=1
		end,
		getInitCondCode = function(self, solver)
			return [[
	const real t = 0.;
	real theta = 2. * M_PI / d * (xc.x - t);
	real H = 1. + initCond->A * sin(theta);
	alpha = sqrt(H);
	gamma_ll.xx = H;
	K_ll.xx = -M_PI * initCond->A / d * cos(theta) / alpha;
]]
		end,
	},
	{
		name = 'testbed - gauge wave - diagonal',
		init = function(self, solver)
			EinsteinInitCond.init(self, solver, args)
			
			error"finishme"
		end,
	},
	{
		name = 'testbed - linear wave',
		createInitStruct = function(self, solver)
			EinsteinInitCond.createInitStruct(self, solver)

-- TODO changing the range upon init causes something to freeze up ...
--			solver.mins = vec3d(-.5, -.5, -.5)
--			solver.maxs = vec3d(-.5, -.5, -.5)
			solver:setBoundaryMethods'periodic'
			self:addGuiVars{
				{name='A', value=1e-8},
				{name='d', value=1},
			}
		end,
		getInitCondCode = function(self, solver)
			return [[
	const real t = 0.;
	real theta = 2. * M_PI / d * (xc.x - t);
	real b = initCond->A * sin(theta);
	gamma_ll.yy += b;
	gamma_ll.zz -= b;
	real db_dt = -2. * M_PI * initCond->A / d * cos(theta);
	K_ll.yy = .5 * db_dt;
	K_ll.zz = -.5 * db_dt;
]]
		end,
	},
	{
		name = 'testbed - linear wave - diagonal',
		init = function(self, solver)
			EinsteinInitCond.init(self, solver, args)
			
			error"finishme"
		end,
	},
	{
		name = 'testbed - Gowdy',
		init = function(self, solver)
			EinsteinInitCond.init(self, solver, args)
			
			error'finishme'
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
			return initEinstein{
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
	return class(EinsteinInitCond, cl)
end)

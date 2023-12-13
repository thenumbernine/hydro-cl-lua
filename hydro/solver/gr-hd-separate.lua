--[[
combining a GR solver (probably BSSNOK-finite-difference + backwards-Euler integrator)
and a HD solver (Roe)
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local ig = require 'imgui'
local vec3d = require 'vec-ffi.vec3d'
local vec3sz = require 'vec-ffi.vec3sz'
local template = require 'template'
local clnumber = require 'cl.obj.number'

--[[
TODO take the composite behavior stuff that this and TwoFluidEMHDBehavior have in common
and make a separate parent out of them ... CompositeSolver or something

also TODO ... think of a better way to provide separate arguments to each sub-solver's initialization
--]]
local GRHDSeparateSolver = class()

GRHDSeparateSolver.name = 'GR+HD'

function GRHDSeparateSolver:init(args)
	local solverargs = args or {}	
	local einsteinargs = solverargs.einstein or {}
	self.app = assert(args.app)

	local GRSolver = require 'hydro.solver.z4c-fd':subclass()
	function GRSolver:init(args)
		GRSolver.super.init(self, table(args, {
			initCond = einsteinargs.initCond or 'Minkowski',
			integrator = einsteinargs.integrator or 'backward Euler',
		}))
		self.name = 'GR '..self.name
	end
	local gr = GRSolver(args)
	self.gr = gr

	local HydroSolver = require 'hydro.solver.grhd-roe':subclass()
	-- TODO :createCodePrefix() has been replaced with :initCodeModules()
	function HydroSolver:createCodePrefix()
		HydroSolver.super.createCodePrefix(self)
		
		self.codePrefix = table{
			self.modules:getCodeAndHeader(self.sharedModulesEnabled:keys():unpack()),
			gr.eqn:getTypeCode(),
			
			-- this is for calc_exp_neg4phi
			gr.eqn:getCommonFuncCode(),
		}:concat'\n'
	end
	function HydroSolver:init(args)
		args = table(args, {
			-- TODO make initConds objects that accept parameters
			--  and make a spherical object at an arbitrary location
			--initCond = 
		})
		HydroSolver.super.init(self, args)
		self.name = 'HD '..self.name
	end
	function HydroSolver:getADMArgs()
		return template([[,
	global <?=gr.eqn.symbols.cons_t?> const * const grUBuf]], {gr=gr})
	end
	-- args is volatile
	function HydroSolver:getADMVarCode(args)
		args = args or {}
		args.suffix = args.suffix or ''
		args.index = args.index or ('index'..args.suffix)
		args.alpha = args.alpha or ('alpha'..args.suffix)
		args.beta = args.beta or ('beta'..args.suffix)
		args.gamma = args.gamma or ('gamma'..args.suffix)
		args.U = 'grU'..args.suffix
		return template([[
	global <?=gr.eqn.symbols.cons_t?> const * const <?=args.U?> = grUBuf + <?=args.index?>;
	real <?=args.alpha?> = <?=args.U?>->alpha;
	real3 <?=args.beta?> = <?=args.U?>->beta_u;
	//real3s3 gammaHat_ll = coord_g_ll(x); // with x I get some redefinitions, without it I get some undefined x's...
	real3s3 gammaHat_ll = coord_g_ll(cellBuf[index].pos);
	real3s3 gammaBar_ll = real3s3_add(gammaHat_ll, <?=args.U?>->epsilon_ll);
	real exp_4phi = 1. / <?=calc_exp_neg4phi?>(<?=args.U?>);
	real3s3 <?=args.gamma?> = real3s3_real_mul(gammaBar_ll, exp_4phi);
]], {gr=gr, args=args})
	end
	
	HydroSolver.DisplayVar_U = HydroSolver.DisplayVar_U:subclass()
	function HydroSolver.DisplayVar_U:setArgs(kernel)
		HydroSolver.DisplayVar_U.super.setArgs(self, kernel)
		kernel:setArg(4, gr.UBuf)
	end
	
	function HydroSolver:getUBufDisplayVarsArgs()
		local args = HydroSolver.super.getUBufDisplayVarsArgs(self)
		args.extraArgs = args.extraArgs or {}
		table.insert(args.extraArgs, 'global '..gr.eqn.symbols.cons_t..' const * const grUBuf')
		args.extraArgNames = args.extraArgNames or {}
		table.insert(args.extraArgNames, 'grUBuf')
		return args
	end
	function HydroSolver:refreshInitStateProgram()
		HydroSolver.super.refreshInitStateProgram(self)
		self.applyInitCondKernelObj.obj:setArg(2, gr.UBuf)
	end
	function HydroSolver:refreshSolverProgram()
		HydroSolver.super.refreshSolverProgram(self)
		
	-- now all of hydro's kernels need to be given the extra ADM arg
io.stderr:write'WARNING!!! make sure gr.UBuf is initialized first!\n'
		self.calcDTKernelObj.obj:setArg(3, gr.UBuf)
		self.addSourceKernelObj.obj:setArg(2, gr.UBuf)
		self.updatePrimsKernelObj.obj:setArg(2, gr.UBuf)
	end

	local hydro = HydroSolver(args)
	self.hydro = hydro

	self.solvers = table{
		self.gr, 
		self.hydro,
	}
	
	self.displayVars = table():append(self.solvers:mapi(function(solver) return solver.displayVars end):unpack())
	
	-- make names unique so that stupid 1D var name-matching code doesn't complain
	self.solverForDisplayVars = table()
	for i,solver in ipairs(self.solvers) do
		for _,var in ipairs(solver.displayVars) do
			self.solverForDisplayVars[var] = solver
			--var.name = solver.name:gsub('[%s]', '_')..'_'..var.name
			var.name = i..'_'..var.name
		end
	end

	-- make a lookup of all vars
	self.displayVarForName = self.displayVars:mapi(function(var)
		return var, var.name
	end)

	self.color = vec3d(math.random(), math.random(), math.random()):normalize()

	self.numGhost = self.hydro.numGhost
	self.dim = self.hydro.dim
	self.gridSize = vec3sz(self.hydro.gridSize)
	self.localSize = vec3sz(self.hydro.localSize:unpack())
	self.globalSize = vec3sz(self.hydro.globalSize:unpack())
	self.sizeWithoutBorder = vec3sz(self.hydro.sizeWithoutBorder)
	self.mins = vec3(self.hydro.mins:unpack())
	self.maxs = vec3(self.hydro.maxs:unpack())

	self.coord = self.hydro.coord
	self.eqn = {
		numStates = self.solvers:mapi(function(solver) return solver.eqn.numStates end):sum(),
		numIntStates = self.solvers:mapi(function(solver) return solver.eqn.numIntStates end):sum(),
		numWaves = self.solvers:mapi(function(solver) return solver.eqn.numWaves end):sum(),
	}

	-- call this after we've assigned 'self' all its fields
	self:replaceSourceKernels()

	self:createCalcStressEnergyKernel()

	self.t = 0
end

function GRHDSeparateSolver:getConsLRTypeCode() return '' end

function GRHDSeparateSolver:createCalcStressEnergyKernel()
	-- build self.codePrefix
	-- TODO FIXME not working
	require 'hydro.solver.gridsolver'.createCodePrefix(self)

	local lines = table{
		self.modules:getCodeAndHeader(self.sharedModulesEnabled:keys():unpack()),
		self.gr.eqn:getTypeCode(),
		self.hydro.eqn:getTypeCode(),
		template([[
kernel void calcStressEnergy(
	global <?=hydro.eqn.symbols.cons_t?>* hydroUBuf,
	global <?=gr.eqn.symbols.cons_t?>* grUBuf
) {
	<?=SETBOUNDS?>(0,0);
	
	global <?=hydro.eqn.symbols.cons_only_t?> * const hydroU = &hydroUBuf[index].cons;
	global <?=hydro.eqn.symbols.prim_t?> * const hydroPrim = &hydroUBuf[index].prim;
	global <?=gr.eqn.symbols.cons_t?> const * const grU = grUBuf + index;

	//as long as t is a separate coordinate ...
	real alpha = grU->alpha;
	real3 beta_u = grU->beta_u;

/*
grU->rho = n^a * n^b * T_ab
grU->S_u = -gamma^ij n^a T_aj
grU->S_ll = gamma_i^c gamma_j^d T_cd

is projection gamma^ij = g^ij + n^i n^j equal to the spatial metric inverse gamma^ij ?
(g^ij + n^i n^j)
= (gamma^ij - beta^i beta^j / alpha^2) + beta^i beta^j / alpha^2
= gamma^ij
yep.


ok so what's our T_ab going to be ...
relativistic fluid: T_ab = rho_mass h u_a u_b + P g_ab 
so rho_gr = n^a n^b T_ab
	= n^a n^b (rho_mass h u_a u_b + P g_ab) 
n^a = (1/alpha, -beta^i/alpha), n_a = (-alpha, 0)
n^a n^b g_ab = n^a n_a = -1 
n_a u^a = -alpha u^t = -W
rho_gr = (n^a u_a)^2 rho_mass h - P
so rho_gr = W^2 rho_mass h - P
*/
	real W = hydroU->prim.W;
	real rhoMass = hydroU->prim.rho;
	real eInt = hydroU->prim.eInt;
	real P = calc_P(rho, eInt);
	real h = calc_h(rhoMass, P, eInt);
	real rhoGR = W * W * rhoMass * h - P;
	grU->rho = rhoGR;

/*
S^i = -gamma^ij n^a T_aj
= -gamma^ij n^a (rho_mass h u_a u_j + P g_aj)
= W rho_mass h u^i
*/
	real3 u_l = real3_add(
		real3_real_mul(hydroU->prim.v, W),
		real3_real_mul(beta_u, -W / alpha));
	real3 u_u = real3s3_mul(gamma_uu, u_l);
	grU->S_u = real3_real_mul(u_u, W * rhoMass * h);

/*
S_ij = gamma_i^a gamma_j^b T_cd
= (delta_i^a + n_i n^a) (delta_j^b + n_j n^b) (rho_mass h u_a u_b + P g_ab)
= rho_mass h u_i u_j 
	+ P g_ij 
	+ P n_i n_j
	- W rho_mass h u_i n_j
	- W rho_mass h n_i u_j 
	+ n_i n_j rho_mass h W^2 
= rho_mass h u_i u_j + P gamma_ij 
= rho_mass h u_i u_j + P gamma_ij 
*/
	grU->S_ll = real3s3_add(
		real3s3_real_mul(real3_outer(u_l), rhoMass * h),
		real3s3_real_mul(gamma_ll, P));
}
]], 	{
			hydro = self.hydro,
			gr = self.gr,
		}),
	}
	local code = lines:concat'\n'
	self.calcStressEnergyProgramObj = self.gr.Program{name='calcStressEnergy', code=code}
	self.calcStressEnergyProgramObj:compile()
	self.calcStressEnergyKernelObj = self.calcStressEnergyProgramObj:kernel('calcStressEnergy', self.hydro.UBuf, self.gr.UBuf)
end

function GRHDSeparateSolver:replaceSourceKernels()

--[=[ instead of copying vars from nr to grhd, I've integrated the nr code directly to the grhd solver
	
	local lines = table{
		self.modules:getCodeAndHeader(self.sharedModulesEnabled:keys():unpack()),
		self.gr.eqn:getTypeCode(),
		self.hydro.eqn:getTypeCode(),
		template([[
kernel void copyMetricFromGRToHydro(
	// I can't remember which of these two uses it -- maybe both? 
	//TODO either just put it in one place ...
	// or better yet, just pass gr.UBuf into all the grhd functions?
	global <?=hydro.eqn.symbols.cons_t?>* hydroUBuf,
	global <?=gr.eqn.symbols.cons_t?> const * const grUBuf
) {
	<?=SETBOUNDS?>(0,0);
	global <?=hydro.eqn.symbols.cons_only_t?>* hydroU = &hydroUBuf[index].cons;
	global <?=hydro.eqn.symbols.prim_t?>* hydroPrim = &hydroUBuf[index].prim;
	global <?=gr.eqn.symbols.cons_t?> const * const grU = grUBuf + index;
	hydroU->alpha = hydroPrim->alpha = grU->alpha;
	hydroU->beta = hydroPrim->beta = grU->beta_u;

	real exp_4phi = exp(4. * grU->phi);
	real3s3 gamma_ll = real3s3_real_mul(grU->gammaTilde_ll, exp_4phi);
	hydroU->gamma = hydroPrim->gamma = gamma_ll;
}
]], 	{
			hydro = self.hydro,
			gr = self.gr,
		}),
	}
	local code = lines:concat'\n'
	self.copyMetricFromGRToHydroProgramObj = self.gr.Program{name='copyMetricFromGRToHydro', code=code}
	self.copyMetricFromGRToHydroProgramObj:compile()
	self.copyMetricFromGRToHydroKernelObj = self.copyMetricFromGRToHydroProgramObj:kernel('copyMetricFromGRToHydro', self.hydro.UBuf, self.gr.UBuf)
-- instead of copying vars from nr to grhd, I've integrated the nr code directly to the grhd solver ]=]
end

function GRHDSeparateSolver:callAll(name, ...)
	local args = setmetatable({...}, table)
	args.n = select('#',...)
	return self.solvers:mapi(function(solver)
		return solver[name](solver, args:unpack(1, args.n))
	end):unpack()
end

function GRHDSeparateSolver:createEqn()
	self:callAll'createEqn'
end

function GRHDSeparateSolver:resetState()
	self:callAll'resetState'
	self.t = self.hydro.t
end

function GRHDSeparateSolver:boundary()
	self:callAll'boundary'
end

function GRHDSeparateSolver:initDraw()
	self:callAll'initDraw'
end

function GRHDSeparateSolver:calcDT()
	return math.min(self:callAll'calcDT')
end

function GRHDSeparateSolver:step(dt)
	-- copy spacetime across for the GRHD
--[=[ instead of copying vars from nr to grhd, I've integrated the nr code directly to the grhd solver
	self.copyMetricFromGRToHydroKernelObj()
-- instead of copying vars from nr to grhd, I've integrated the nr code directly to the grhd solver ]=]
-- [=[ TODO eventually inline the stress-energy calcs.  until then ...
	self.calcStressEnergyKernelObj()
--]=]	
	self:callAll('step', dt)
end

-- same as Solver.update
-- note this means sub-solvers' update() will be skipped
-- so best to put update stuff in step()
function GRHDSeparateSolver:update()
	local dt = self:calcDT()
	self:step(dt)
	self.t = self.hydro.t
end

function GRHDSeparateSolver:getTex(var) 
	return self.solverForDisplayVars[var].tex
end

function GRHDSeparateSolver:calcDisplayVarToTex(var)
	return self.solverForDisplayVars[var]:calcDisplayVarToTex(var)
end

function GRHDSeparateSolver:calcDisplayVarRange(var)
	return self.solverForDisplayVars[var]:calcDisplayVarRange(var)
end

function GRHDSeparateSolver:updateGUI()
	for i,solver in ipairs(self.solvers) do
		ig.igPushID_Str('subsolver '..i)
		if ig.igCollapsingHeader('sub-solver '..solver.name..':') then
			solver:updateGUI()
		end
		ig.igPopID()
	end
end

return GRHDSeparateSolver

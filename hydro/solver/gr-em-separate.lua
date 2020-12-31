--[[
GR + EM
copied from GR + HD
I should make a parent class/template
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local ig = require 'ffi.imgui'
local vec3d = require 'vec-ffi.vec3d'
local vec3sz = require 'vec-ffi.vec3sz'
local template = require 'template'
local clnumber = require 'cl.obj.number'

local GREMSeparateSolver = class()

GREMSeparateSolver.name = 'GR+EM'

function GREMSeparateSolver:init(args)
	self.app = assert(args.app)

	local GRSolver = class(require 'hydro.solver.bssnok-fd')
	function GRSolver:init(args)
		GRSolver.super.init(self, table(args, {
			--initCond = 'black hole - isotropic',
			-- flat-space initial state (how come there is dynamics with this?)
			integrator = 'backward Euler',
		}))
		self.name = 'GR '..self.name
	end
	local gr = GRSolver(args)
	self.gr = gr

	local GRMaxwellSolver = class(require 'hydro.solver.gr-maxwell-roe')
	-- TODO :createCodePrefix() has been replaced with :initCodeModules()
	function GRMaxwellSolver:createCodePrefix()
		GRMaxwellSolver.super.createCodePrefix(self)
		
		self.codePrefix = table{
			self.modules:getCodeAndHeader(self.sharedModulesEnabled:keys():unpack()),
			gr.eqn:getTypeCode(),
			
			-- this is for gr's calc_exp_neg4phi, which em will need 
			gr.eqn:getCommonFuncCode(),
		}:concat'\n'
	end
	function GRMaxwellSolver:init(args)
		args = table(args, {
			eqn = 'gr-maxwell',
		})
		GRMaxwellSolver.super.init(self, args)
		self.name = 'EM '..self.name
	end
	function GRMaxwellSolver:getADMArgs()
		return template([[,
	global <?=gr.eqn.symbols.cons_t?> const * const grUBuf]], {gr=gr})
	end
	-- args is volatile
	function GRMaxwellSolver:getADMVarCode(args)
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
	sym3 <?=args.gamma?> = sym3_real_mul(<?=args.U?>->gammaTilde_ll, 1. / <?=calc_exp_neg4phi?>(<?=args.U?>));
]], {gr=gr, args=args})
	end

	GRMaxwellSolver.DisplayVar_U = class(GRMaxwellSolver.DisplayVar_U)
	function GRMaxwellSolver.DisplayVar_U:setArgs(kernel)
		GRMaxwellSolver.DisplayVar_U.super.setArgs(self, kernel)
		kernel:setArg(4, gr.UBuf)
	end
	
	function GRMaxwellSolver:getUBufDisplayVarsArgs()
		local args = GRMaxwellSolver.super.getUBufDisplayVarsArgs(self)
		args.extraArgs = args.extraArgs or {}
		table.insert(args.extraArgs, 'global '..gr.eqn.symbols.cons_t..' const * const grUBuf')
		args.extraArgNames = args.extraArgNames or {}
		table.insert(args.extraArgNames, 'grUBuf')
		return args
	end
	function GRMaxwellSolver:refreshInitStateProgram()
		GRMaxwellSolver.super.refreshInitStateProgram(self)
		self.applyInitCondKernelObj.obj:setArg(2, gr.UBuf)
	end
	function GRMaxwellSolver:refreshSolverProgram()
		GRMaxwellSolver.super.refreshSolverProgram(self)
		
		-- now all of em's kernels need to be given the extra ADM arg
		-- TODO replace all these kernels with einstein-maxwell equivalents
		-- in fact, a better option might be like grhd is for euler ...
		-- ... just create a einstein-maxwell file that is missing alpha,beta,gamma variables
		-- and have the GR solver provide that code
io.stderr:write'WARNING!!! make sure gr.UBuf is initialized first!\n'
		self.calcDTKernelObj.obj:setArg(3, gr.UBuf)
		self.addSourceKernelObj.obj:setArg(2, gr.UBuf)
	end

	local em = GRMaxwellSolver(args)
	self.em = em

	self.solvers = table{
		self.gr, 
		self.em,
	}
	
	self.displayVars = table():append(self.solvers:map(function(solver) return solver.displayVars end):unpack())
	
	-- make names unique so that stupid 1D var name-matching code doesn't complain
	self.solverForDisplayVars = table()
	for i,solver in ipairs(self.solvers) do
		for _,var in ipairs(solver.displayVars) do
			self.solverForDisplayVars[var] = solver
			--var.name = solver.name:gsub('[%s]', '_')..'_'..var.name
			var.name = i..'_'..var.name
		end
	end

	self.color = vec3d(math.random(), math.random(), math.random()):normalize()

	self.numGhost = self.em.numGhost
	self.dim = self.em.dim
	self.gridSize = vec3sz(self.em.gridSize)
	self.localSize = vec3sz(self.em.localSize:unpack())
	self.globalSize = vec3sz(self.em.globalSize:unpack())
	self.sizeWithoutBorder = vec3sz(self.em.sizeWithoutBorder)
	self.mins = vec3d(self.em.mins:unpack())
	self.maxs = vec3d(self.em.maxs:unpack())

	self.coord = self.em.coord
	self.eqn = {
		numStates = self.solvers:map(function(solver) return solver.eqn.numStates end):sum(),
		numIntStates = self.solvers:map(function(solver) return solver.eqn.numIntStates end):sum(),
		numWaves = self.solvers:map(function(solver) return solver.eqn.numWaves end):sum(),
	}

	-- call this after we've assigned 'self' all its fields
	self:replaceSourceKernels()

	self.t = 0
end

function GREMSeparateSolver:getConsLRTypeCode() return '' end

function GREMSeparateSolver:replaceSourceKernels()

	local lines = table{
		self.modules:getCodeAndHeader(self.sharedModulesEnabled:keys():unpack()),
		self.gr.eqn:getTypeCode(),
		-- TODO :getCommonFuncCode() has been replaced with :initCodeModules()
		self.gr.eqn:getCommonFuncCode(),
		self.em.eqn:getTypeCode(),
		template([[
kernel void computeGRStressEnergy(
	global <?=gr.eqn.symbols.cons_t?> * const grUBuf,
	global <?=em.eqn.symbols.cons_t?> const * const emUBuf
) {
	<?=SETBOUNDS?>(solver->numGhost, solver->numGhost);
	
	//populate grUBuf->rho, S_u, S_ll
	
	global <?=gr.eqn.symbols.cons_t?> * const grU = grUBuf + index;
	global <?=em.eqn.symbols.cons_t?> const * const emU = emUBuf + index;

	real exp_neg4phi = <?=calc_exp_neg4phi?>(grU);
	sym3 gamma_ll = sym3_real_mul(grU->gammaTilde_ll, 1. / exp_neg4phi);
	sym3 gamma_uu = sym3_real_mul(grU->gammaTilde_uu, exp_neg4phi);
	//g^ij = gamma^ij - beta^i beta^j / alpha^2
	sym3 g_uu = sym3_sub(
		gamma_uu, 
		real3_outer(
			real3_real_mul(
				grU->beta_u, 
				1. / grU->alpha)));
	
	real3 B_u = emU->B;
	real3 B_l = sym3_real3_mul(gamma_ll, B_u);
	real3 E_u = real3_real_mul(emU->epsE, 1. / emU->eps);
	real3 E_l = sym3_real3_mul(gamma_ll, E_u);

	//2009 Alcubierre eqn 2.75-2.77
	real ESqBSq = real3_weightedLenSq(E_u, gamma_ll) + real3_weightedLenSq(B_u, gamma_ll);
	grU->rho = ESqBSq * (1. / (8. * M_PI));

	//S_j = [jlm] E^l B^m
	//g^ab J_b = g^it J_t + g^ij J_j = g^ij S_j / (4 pi)
	//now does this raise by g^ij and not gamma^ij?  the eqn looks like so...
	real3 S_l = real3_cross(E_u, B_u);
	grU->S_u = real3_real_mul(sym3_real3_mul(g_uu, S_l), 1. / (4. * M_PI));

	//E_i = g_ia E^a = g_ij E^j (since E^t = 0) = gamma_ij E^j (by ADM metric def)
	//S_ij = (gamma_ij (E^2 + B^2) - 2 (E_i E_j + B_i B_j)) / (8 pi)
	grU->S_ll = sym3_sub(
		sym3_real_mul(gamma_ll, ESqBSq / (8. * M_PI)),
		sym3_real_mul(
			sym3_add(
				real3_outer(E_l),
				real3_outer(B_l)
			), 1. / (4. * M_PI)
		)
	);
}
]], 		{
			em = self.em,
			gr = self.gr,
		}),
	}
	local code = lines:concat'\n'
	self.computeGRStressEnergyProgramObj = self.gr.Program{name='computeGRStressEnergy', code=code}
	self.computeGRStressEnergyProgramObj:compile()
	self.computeGRStressEnergyKernelObj = self.computeGRStressEnergyProgramObj:kernel('computeGRStressEnergy', self.gr.UBuf, self.em.UBuf)
end

function GREMSeparateSolver:callAll(name, ...)
	local args = setmetatable({...}, table)
	args.n = select('#',...)
	return self.solvers:map(function(solver)
		return solver[name](solver, args:unpack(1, args.n))
	end):unpack()
end

function GREMSeparateSolver:createEqn()
	self:callAll'createEqn'
end

function GREMSeparateSolver:resetState()
	self:callAll'resetState'
	self.t = self.em.t
end

function GREMSeparateSolver:boundary()
	self:callAll'boundary'
end

local function passthru(f, ...)
	f(...)
	return ...
end

function GREMSeparateSolver:calcDT()
	return math.min(self:callAll'calcDT')
end

function GREMSeparateSolver:step(dt)
	-- this computes rho, S_u, S_ll for the GR solver
	--self.computeGRStressEnergyKernelObj()
	
	self:callAll('step', dt)
end

-- same as Solver.update
-- note this means sub-solvers' update() will be skipped
-- so best to put update stuff in step()
function GREMSeparateSolver:update()
	local dt = self:calcDT()
	self:step(dt)
	self.t = self.em.t
end

function GREMSeparateSolver:getTex(var) 
	return self.solverForDisplayVars[var].tex
end

function GREMSeparateSolver:calcDisplayVarToTex(var)
	return self.solverForDisplayVars[var]:calcDisplayVarToTex(var)
end

function GREMSeparateSolver:calcDisplayVarRange(var)
	return self.solverForDisplayVars[var]:calcDisplayVarRange(var)
end

function GREMSeparateSolver:updateGUI()
	for i,solver in ipairs(self.solvers) do
		ig.igPushIDStr('subsolver '..i)
		if ig.igCollapsingHeader('sub-solver '..solver.name..':') then
			solver:updateGUI()
		end
		ig.igPopID()
	end
end

return GREMSeparateSolver

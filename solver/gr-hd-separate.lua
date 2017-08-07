--[[
combining a GR solver (probably BSSNOK-finite-difference + backwards-Euler integrator)
and a HD solver (Roe)
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local ig = require 'ffi.imgui'
local vec3sz = require 'ffi.vec.vec3sz'
local template = require 'template'
local CLProgram = require 'cl.program'
local clnumber = require 'cl.obj.number'

--[[
TODO take the composite behavior stuff that this and TwoFluidEMHDBehavior have in common
and make a separate parent out of them ... CompositeSolver or something

also TODO ... think of a better way to provide separate arguments to each sub-solver's initialization
--]]
local GRHDSeparateSolver = class()

GRHDSeparateSolver.name = 'GR+HD'

function GRHDSeparateSolver:init(args)
	self.app = assert(args.app)

	local GRSolver = class(require 'solver.bssnok-fd')
	function GRSolver:init(args)
		GRSolver.super.init(self, table(args, {
			initState = 'black hole - isotropic',
			integrator = 'backward Euler',
		}))
		self.name = 'GR '..self.name
	end
	local gr = GRSolver(args)
	self.gr = gr

	local HydroSolver = class(require 'solver.grhd-roe')
	function HydroSolver:createCodePrefix()
		HydroSolver.super.createCodePrefix(self)
		self.codePrefix = table{
			self.codePrefix,
			gr.eqn:getTypeCode(),
		}:concat'\n'
	end
	function HydroSolver:init(args)
		args = table(args, {
			-- TODO make initStates objects that accept parameters
			--  and make a spherical object at an arbitrary location
			--initState = 
		})
		HydroSolver.super.init(self, args)
		self.name = 'HD '..self.name
	end
	function HydroSolver:getADMArgs()
		return template([[,
	const global <?=gr.eqn.cons_t?>* grUBuf]], {gr=gr})
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
	const global <?=gr.eqn.cons_t?>* <?=args.U?> = grUBuf + <?=args.index?>;
	real <?=args.alpha?> = <?=args.U?>->alpha;
	real3 <?=args.beta?> = <?=args.U?>->beta_u;
	sym3 <?=args.gamma?> = sym3_scale(<?=args.U?>->gammaBar_ll, exp(4. * <?=args.U?>->phi));
]], {gr=gr, args=args})
	end
	
	HydroSolver.ConvertToTex_U = class(HydroSolver.ConvertToTex_U)
	function HydroSolver.ConvertToTex_U:setArgs(kernel, var)
		HydroSolver.ConvertToTex_U.super.setArgs(self, kernel, var)
		kernel:setArg(3, gr.UBuf)
	end
	function HydroSolver:getAddConvertToTexUBufArgs()
		local args = HydroSolver.super.getAddConvertToTexUBufArgs(self)
		table.insert(args.extraArgs, 'const global '..gr.eqn.cons_t..'* grUBuf')
		return args
	end
	function HydroSolver:refreshInitStateProgram()
		HydroSolver.super.refreshInitStateProgram(self)
		self.initStateKernel:setArg(2, gr.UBuf)
	end
	function HydroSolver:refreshSolverProgram()
		HydroSolver.super.refreshSolverProgram(self)
		
		self.calcDTKernel:setArg(2, gr.UBuf)
		self.calcEigenBasisKernel:setArg(3, gr.UBuf)
	end

	local hydro = HydroSolver(args)
	self.hydro = hydro

	-- now all of hydro's kernels need to be given the extra ADM arg
io.stderr:write'WARNING!!! make sure gr.UBuf is initialized first!\n'
	hydro.initStateKernel:setArg(2, gr.UBuf)
	hydro.calcDTKernel:setArg(2, gr.UBuf)
	hydro.calcEigenBasisKernel:setArg(3, gr.UBuf)
	hydro.updatePrimsKernel:setArg(2, gr.UBuf)

	self.solvers = table{
		self.gr, 
		self.hydro,
	}
	
	self.displayVars = table():append(self.solvers:map(function(solver) return solver.displayVars end):unpack())
	
	-- make names unique so that stupid 1D var name-matching code doesn't complain
	self.solverForDisplayVars = table()
	for _,solver in ipairs(self.solvers) do
		for _,var in ipairs(solver.displayVars) do
			self.solverForDisplayVars[var] = solver
			var.name = solver.name:gsub('[%s]', '_')..'_'..var.name
		end
	end

	self.color = vec3(math.random(), math.random(), math.random()):normalize()

	self.numGhost = self.hydro.numGhost
	self.dim = self.hydro.dim
	self.gridSize = vec3sz(self.hydro.gridSize)
	self.localSize = vec3sz(self.hydro.localSize:unpack())
	self.globalSize = vec3sz(self.hydro.globalSize:unpack())
	self.sizeWithoutBorder = vec3sz(self.hydro.sizeWithoutBorder)
	self.mins = vec3(self.hydro.mins:unpack())
	self.maxs = vec3(self.hydro.maxs:unpack())

	self.dxs = vec3(self.hydro.dxs:unpack())
	self.geometry = self.hydro.geometry
	self.eqn = {
		numStates = self.solvers:map(function(solver) return solver.eqn.numStates end):sum(),
		numWaves = self.solvers:map(function(solver) return solver.eqn.numWaves end):sum(),
		getEigenTypeCode = function() end,
		getCodePrefix = function() end,
	}

	-- call this after we've assigned 'self' all its fields
	self:replaceSourceKernels()

	self.t = 0
end

function GRHDSeparateSolver:getConsLRTypeCode() return '' end

function GRHDSeparateSolver:replaceSourceKernels()

--[=[ instead of copying vars from nr to grhd, I've integrated the nr code directly to the grhd solver
	
	-- build self.codePrefix
	require 'solver.solver'.createCodePrefix(self)
	
	local lines = table{
		self.app.env.code,
		self.codePrefix,
		self.gr.eqn:getTypeCode(),
		self.hydro.eqn:getTypeCode(),
		template([[
kernel void copyMetricFromGRToHydro(
// I can't remember which of these two uses it -- maybe both? 
//TODO either just put it in one place ...
// or better yet, just pass gr.UBuf into all the grhd functions?
global <?=hydro.eqn.cons_t?>* hydroUBuf,
global <?=hydro.eqn.prim_t?>* hydroPrimBuf,
const global <?=gr.eqn.cons_t?>* grUBuf
) {
SETBOUNDS(0,0);
global <?=hydro.eqn.cons_t?>* hydroU = hydroUBuf + index;
global <?=hydro.eqn.prim_t?>* hydroPrim = hydroPrimBuf + index;
const global <?=gr.eqn.cons_t?>* grU = grUBuf + index;
hydroU->alpha = hydroPrim->alpha = grU->alpha;
hydroU->beta = hydroPrim->beta = grU->beta_u;

real exp_4phi = exp(4. * grU->phi);
sym3 gamma_ll = sym3_scale(grU->gammaBar_ll, exp_4phi);
hydroU->gamma = hydroPrim->gamma = gamma_ll;
}
]], 		{
			hydro = self.hydro,
			gr = self.gr,
		}),
	}
	local code = lines:concat'\n'
	self.copyMetricFrmoGRToHydroProgram = CLProgram{context=self.app.ctx, devices={self.app.device}, code=code}
	self.copyMetricFrmoGRToHydroKernel = self.copyMetricFrmoGRToHydroProgram:kernel('copyMetricFromGRToHydro', self.hydro.UBuf, self.hydro.primBuf, self.gr.UBuf)
-- instead of copying vars from nr to grhd, I've integrated the nr code directly to the grhd solver ]=]
end

function GRHDSeparateSolver:callAll(name, ...)
	local args = setmetatable({...}, table)
	args.n = select('#',...)
	return self.solvers:map(function(solver)
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

function GRHDSeparateSolver:calcDT()
	return math.min(self:callAll'calcDT')
end

function GRHDSeparateSolver:step(dt)
	-- copy spacetime across for the GRHD
--[=[ instead of copying vars from nr to grhd, I've integrated the nr code directly to the grhd solver
	self.app.cmds:enqueueNDRangeKernel{kernel=self.copyMetricFrmoGRToHydroKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
-- instead of copying vars from nr to grhd, I've integrated the nr code directly to the grhd solver ]=]
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
		ig.igPushIdStr('subsolver '..i)
		if ig.igCollapsingHeader('sub-solver '..solver.name..':') then
			solver:updateGUI()
		end
		ig.igPopId()
	end
end

return GRHDSeparateSolver

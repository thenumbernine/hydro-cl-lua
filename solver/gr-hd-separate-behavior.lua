--[[
combining a GR solver (probably BSSNOK-finite-difference + backwards-Euler integrator)
and a HD solver (Roe)

this is no longer a 'behavior' but now fully a 'solver'
though only an abstract version of one, hiding the two underlying solvers
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

also TODO ... think of a better way to provide separate arguments to each sub-solver

right now 'parent' is just the Euler solver type
I'm only using BSSNOK-finite-difference for the spacetime solver
--]]
local function GRHDBehavior()
	local templateClass = class()

	templateClass.name = 'GR+HD'

	function templateClass:init(args)
		self.app = assert(args.app)

		local HydroSolver = class(require 'solver.grhd-roe')
		function HydroSolver:init(args)
			HydroSolver.super.init(self, args)
			self.name = 'HD '..self.name
		end
		self.hydro = HydroSolver(args)

		local GRSolver = class(require 'solver.bssnok-fd')
		function GRSolver:init(args)
			GRSolver.super.init(self, table(args, {
				initState = 'black hole - isotropic',
				integrator = 'backward Euler',
			}))
			self.name = 'GR '..self.name
		end
		self.gr = GRSolver(args)

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

	function templateClass:getConsLRTypeCode() return '' end

	function templateClass:replaceSourceKernels()
		
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
		self.copyMetricFrmoGRToHydroKernel = self.copyMetricFrmoGRToHydroProgram:kernel(
			'copyMetricFromGRToHydro',
			self.hydro.UBuf,
			self.hydro.primBuf,
			self.gr.UBuf)
	end
	
	function templateClass:callAll(name, ...)
		local args = setmetatable({...}, table)
		args.n = select('#',...)
		return self.solvers:map(function(solver)
			return solver[name](solver, args:unpack(1, args.n))
		end):unpack()
	end
	
	function templateClass:createEqn()
		self:callAll'createEqn'
	end
	
	function templateClass:resetState()
		self:callAll'resetState'
		self.t = self.hydro.t
	end
	
	function templateClass:boundary()
		self:callAll'boundary'
	end
	
	function templateClass:calcDT()
		return math.min(self:callAll'calcDT')
	end
	
	function templateClass:step(dt)
		self:callAll('step', dt)
		self.app.cmds:enqueueNDRangeKernel{kernel=self.copyMetricFrmoGRToHydroKernel, dim=self.dim, globalSize=self.globalSize:ptr(), localSize=self.localSize:ptr()}
	end
	
	-- same as Solver.update
	-- note this means sub-solvers' update() will be skipped
	-- so best to put update stuff in step()
	function templateClass:update()
		local dt = self:calcDT()
		self:step(dt)
		self.t = self.hydro.t
	end

	function templateClass:getTex(var) 
		return self.solverForDisplayVars[var].tex
	end

	function templateClass:calcDisplayVarToTex(var)
		return self.solverForDisplayVars[var]:calcDisplayVarToTex(var)
	end

	function templateClass:calcDisplayVarRange(var)
		return self.solverForDisplayVars[var]:calcDisplayVarRange(var)
	end

	function templateClass:updateGUI()
		for i,solver in ipairs(self.solvers) do
			ig.igPushIdStr('subsolver '..i)
			if ig.igCollapsingHeader('sub-solver '..solver.name..':') then
				solver:updateGUI()
			end
			ig.igPopId()
		end
	end

	return templateClass
end

return GRHDBehavior

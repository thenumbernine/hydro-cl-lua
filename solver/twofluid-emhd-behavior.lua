--[[
using 2014 Abgrall, Kumar "Robust Finite Volume Schemes for Two-Fluid Plasma Equations

mass_ion / mass_electron is 1836 in real life
paper uses 25 
also puts ion pressure at 1/100'th the electron pressure

--]]

local class = require 'ext.class'
local table = require 'ext.table'
local ig = require 'ffi.imgui'
local vec3sz = require 'ffi.vec.vec3sz'
local template = require 'template'

--[[
parent is a solver (Roe, RoeImplicitLinearized)
--]]
local function TwoFluidEMHDBehavior(parent)
	local templateClass = class()

	templateClass.name = 'two-fluid EMHD '..parent.name

	function templateClass:init(args)
		self.app = assert(args.app)

		-- how to specify initial conditions for the both of them?
		-- separate args? maxwellInitConds vs eulerInitConds?
		-- same name in init/euler and init/maxwell?
		-- both?

		local IonSolver = class(require 'solver.euler-behavior'(parent))
		function IonSolver:init(args)
			IonSolver.super.init(self, table(args, {
--				initState = 'two-fluid EMHD soliton ion',
			}))
			self.name = 'ion '..self.name
		end
		self.ion = IonSolver(args)

		local ElectronSolver = class(require 'solver.euler-behavior'(parent))
		function ElectronSolver:init(args)
			ElectronSolver.super.init(self, table(args, {
--				initState = 'two-fluid EMHD soliton electron',
			}))
			self.name = 'electron '..self.name
		end
		self.electron = ElectronSolver(args)

		local MaxwellSolver = class(require 'solver.maxwell-behavior'(parent))
		function MaxwellSolver:init(args)
			MaxwellSolver.super.init(self, table(args, {
--				initState = 'two-fluid EMHD soliton maxwell',
			}))
		end
		self.maxwell = MaxwellSolver(args)
		
		self.solvers = table{self.ion, self.electron, self.maxwell}
		
		self.displayVars = table():append(self.solvers:map(function(solver) return solver.displayVars end):unpack())
		
		-- make names unique so that stupid 1D var name-matching code doesn't complain
		self.solverForDisplayVars = table()
		for _,solver in ipairs(self.solvers) do
			for _,var in ipairs(solver.displayVars) do
				self.solverForDisplayVars[var] = solver
				var.name = solver.name:gsub('[%s]', '_')..'_'..var.name
			end
		end

		self.ion:refreshDisplayProgram()
		self.electron:refreshDisplayProgram()
		self.maxwell:refreshDisplayProgram()

		self.color = vec3(math.random(), math.random(), math.random()):normalize()
		
		self.numGhost = self.ion.numGhost
		self.dim = self.ion.dim
		self.gridSize = vec3sz(self.ion.gridSize)
		self.mins = vec3(self.ion.mins:unpack())
		self.maxs = vec3(self.ion.maxs:unpack())
			
		-- only used by createCodePrefix
		self.dxs = vec3(self.ion.dxs:unpack())
		self.geometry = self.ion.geometry
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

	-- only used by createCodePrefix
	function templateClass:getCoordMapCode() 
		return self.ion:getCoordMapCode() 
	end

	function templateClass:getConsLRTypeCode() return '' end

	local clnumber = require 'clnumber'
	function templateClass:replaceSourceKernels()
		local chargeMassRatio_ion = 1
		local chargeMassRatio_electron = .01
		local eps0 = 1

		require 'solver.roe'.createCodePrefix(self)

		local lines = table{
			self.codePrefix,
			self.electron.eqn:getTypeCode(),
			self.ion.eqn:getTypeCode(),
			self.maxwell.eqn:getTypeCode(),

			'#define chargeMassRatio_ion '..clnumber(chargeMassRatio_ion),
			'#define chargeMassRatio_electron '..clnumber(chargeMassRatio_electron),
			'#define eps0 '..clnumber(eps0),
			[[
<? for _,species in ipairs{'ion', 'electron'} do ?>

kernel void addSource_<?=species?>(
	global <?=euler_cons_t?>* derivBuf,
	const global <?=euler_cons_t?>* UBuf,
	const global <?=maxwell_cons_t?>* maxwellUBuf
) {
	SETBOUNDS(2,2);
	global <?=euler_cons_t?>* deriv = derivBuf + index;
	const global <?=euler_cons_t?>* U = UBuf + index;
	const global <?=maxwell_cons_t?>* maxwellU = maxwellUBuf + index;
	deriv->m.x += chargeMassRatio_<?=species?> * (maxwellU->epsE.x / eps0 + U->m.y * maxwellU->B.z - U->m.z * maxwellU->B.y);
	deriv->m.y += chargeMassRatio_<?=species?> * (maxwellU->epsE.y / eps0 + U->m.z * maxwellU->B.x - U->m.x * maxwellU->B.z);
	deriv->m.z += chargeMassRatio_<?=species?> * (maxwellU->epsE.z / eps0 + U->m.x * maxwellU->B.y - U->m.y * maxwellU->B.x);
	deriv->ETotal += chargeMassRatio_<?=species?> * real3_dot(maxwellU->epsE, U->m) / eps0;
}

<? end ?>

kernel void addSource_maxwell(
	global <?=maxwell_cons_t?>* derivBuf,
	const global <?=euler_cons_t?>* ionUBuf,
	const global <?=euler_cons_t?>* electronUBuf
) {
	SETBOUNDS(2,2);
	global <?=maxwell_cons_t?>* deriv = derivBuf + index;
	const global <?=euler_cons_t?>* ionU = ionUBuf + index;
	const global <?=euler_cons_t?>* electronU = electronUBuf + index;
	deriv->epsE = real3_sub(deriv->epsE, real3_scale(ionU->m, chargeMassRatio_ion));
	deriv->epsE = real3_sub(deriv->epsE, real3_scale(electronU->m, chargeMassRatio_electron));
}

]]
		}
		local code = lines:concat'\n'
		code = template(code, {
			solver = self,
			euler_cons_t = self.electron.eqn.cons_t,
			maxwell_cons_t = self.maxwell.eqn.cons_t,
		})
		self.addSourceProgram = require 'cl.program'{context=self.app.ctx, devices={self.app.device}, code=code}
		
		self.ion.addSourceKernel = self.addSourceProgram:kernel'addSource_ion'
		self.ion.addSourceKernel:setArg(1, self.ion.UBuf)
		self.ion.addSourceKernel:setArg(2, self.maxwell.UBuf)
		self.ion.eqn.useSourceTerm = true
		
		self.electron.addSourceKernel = self.addSourceProgram:kernel'addSource_electron'
		self.electron.addSourceKernel:setArg(1, self.electron.UBuf)
		self.electron.addSourceKernel:setArg(2, self.maxwell.UBuf)
		self.electron.eqn.useSourceTerm = true

		self.maxwell.addSourceKernel = self.addSourceProgram:kernel'addSource_maxwell'
		self.maxwell.addSourceKernel:setArg(1, self.ion.UBuf)
		self.maxwell.addSourceKernel:setArg(2, self.electron.UBuf)
	end

	function templateClass:callAll(name, ...)
		local args = setmetatable({...}, table)
		args.n = select('#',...)
		return self.solvers:map(function(solver)
			return solver[name](solver, args:unpack(1, args.n))
		end):unpack()
	end

	function templateClass:getCoordMapGLSLCode()
		return self.ion:getCoordMapGLSLCode()
	end

	function templateClass:createEqn()
		self:callAll'createEqn'
	end

	function templateClass:resetState()
		self:callAll'resetState'
		self.t = self.ion.t
	end

	function templateClass:boundary()
		self:callAll'boundary'
	end

	function templateClass:calcDT()
		return math.min(self:callAll'calcDT')
	end

	function templateClass:step(dt)
		self:callAll('step', dt)
	end

	-- same as Solver.update
	-- note this means sub-solvers' update() will be skipped
	-- so best to put update stuff in step()
	function templateClass:update()
		local dt = self:calcDT()
		self:step(dt)
		self.t = self.ion.t
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

return TwoFluidEMHDBehavior

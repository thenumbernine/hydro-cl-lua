--[[
using 2014 Abgrall, Kumar "Robust Finite Volume Schemes for Two-Fluid Plasma Equations

mass_ion / mass_electron is 1836 in real life
paper uses 25 
also puts ion pressure at 1/100'th the electron pressure

--]]
local ffi = require 'ffi'
local class = require 'ext.class'
local table = require 'ext.table'
local ig = require 'imgui'
local vec3d = require 'vec-ffi.vec3d'
local vec3sz = require 'vec-ffi.vec3sz'
local template = require 'template'
local clnumber = require 'cl.obj.number'

--[[
parent is a solver... just Roe for now
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

		--[[ specific initialization
		local ionInitState = 'two-fluid EMHD soliton ion'
		local electronInitState = 'two-fluid EMHD soliton electron'
		local emhdInitState = 'two-fluid EMHD soliton maxwell'
		--]]	
		-- [[ use whatever is in args
		local ionInitState
		local electronInitState
		local emhdInitState
		-- TODO scale elec->U.rho and elec->U.P down by 1/ionElectronMassRatio
		--]]

		local IonSolver = parent:subclass()
		IonSolver.eqnName = 'euler'
		function IonSolver:init(args)
			IonSolver.super.init(self, table(args, {initCond = ionInitState}))
			self.name = 'ion '..self.name
		end
		self.ion = IonSolver(args)

		local ElectronSolver = parent:subclass()
		ElectronSolver.eqnName = 'euler'
		function ElectronSolver:init(args)
			ElectronSolver.super.init(self, table(args, {initCond = electronInitState}))
			self.name = 'electron '..self.name
		end
		self.electron = ElectronSolver(args)

		local MaxwellSolver = parent:subclass()
		MaxwellSolver.eqnName = 'glm-maxwell'
		function MaxwellSolver:init(args)
			MaxwellSolver.super.init(self, table(args, {initCond = emhdInitState}))
		end
		self.maxwell = MaxwellSolver(args)

		self.solvers = table{
			self.ion, 
			self.electron,
			self.maxwell,
		}
		
		self.displayVars = table():append(self.solvers:mapi(function(solver) return solver.displayVars end):unpack())
		
		-- make names unique so that stupid 1D var name-matching code doesn't complain
		self.solverForDisplayVars = table()
		for _,solver in ipairs(self.solvers) do
			for _,var in ipairs(solver.displayVars) do
				self.solverForDisplayVars[var] = solver
				var.name = solver.name:gsub('[%s]', '_')..'_'..var.name
			end
		end

		self.color = vec3d(math.random(), math.random(), math.random()):normalize()
		
		self.numGhost = self.ion.numGhost
		self.dim = self.ion.dim
		self.gridSize = vec3sz(self.ion.gridSize)
		self.sizeWithoutBorder = vec3sz(self.ion.sizeWithoutBorder)
		self.mins = vec3d(self.ion.mins:unpack())
		self.maxs = vec3d(self.ion.maxs:unpack())
			
		-- only used by createCodePrefix
		self.coord = self.ion.coord
		self.eqn = {
			numStates = self.solvers:mapi(function(solver) return solver.eqn.numStates end):sum(),
			numIntStates = 0,	-- hmm, nothing should be using this anyways ...
			numWaves = 0,	-- nor this ..
			getTypeCode = function()
				return table{
					'// TwoFluidEMHDSeparateBehavior.eqn:getTypeCode() begin',
					self.electron.eqn:getTypeCode(),
					self.ion.eqn:getTypeCode(),
					self.maxwell.eqn:getTypeCode(),
					'// TwoFluidEMHDSeparateBehavior.eqn:getTypeCode() end',
				}:concat'\n'
			end,
			getEnv = function() return self.maxwell.eqn:getEnv() end,
		}

		setmetatable(self.eqn, require 'hydro.eqn.eqn')
		-- taken from hydro/eqn/twofluid-emhd.lua
		self.eqn:addGuiVars(self.eqn, table{
			--never given, only stated as "speeds for the Maxwell equation"
			-- of course, they're associated with the potentials, so... they could be arbitrary
			{name='divPsiWavespeed', value=1},
			{name='divPhiWavespeed', value=1},
			
			-- gamma = heat capacity ratio
			{name='heatCapacityRatio', value=5/3},

			-- x_0 = reference length
			{name='referenceLength', value=1},

			-- B_0 = reference magnetic field
			{name='referenceMagneticField', value=1},

			-- v_i^T = ion reference thermal velocity
			--{name='ionReferenceThermalVelocity', value=1},
			-- normalized Larmor radius:
			-- lHat = l_r / x_0 = gamma m_i v_i^T / (q_i B_0 x_0)	
			-- the 2014 Abgrall, Kumar model has no gamma, but Wiki says to add gamma = Lorentz boost
			-- therefore:
			-- v_i^T = l_r B_0 q_i / m_i
			-- just use this equation:
			-- v_i^T = l_r r_i B_0

			-- m_i = ion mass
			{name='ionMass', value=1},
		
			-- c = speed of light
			{name='speedOfLight', value=1},
			-- normalized speed of light:
			-- cHat = c / v_i^T
			
			-- l_r = larmor radius
			-- in the experimental section of 2014 Abgrall Kumar this is given fixed values
			{name='ionLarmorRadius', value=.1},
			
			-- m = m_i / m_e
			-- https://en.wikipedia.org/wiki/Proton-to-electron_mass_ratio
			-- m_i / m_e = 1836.15267389
			{name='ionElectronMassRatio', value=100},
			
			-- r_i = q_i / m_i
			{name='ionChargeMassRatio', value=1},
		
			-- lambda_d = ion Debye length
			-- lambda_d = sqrt(epsilon (v_i^T)^2 / n_0 q_i)
			-- the 2014 Abgrall, Kumar model uses vacuum permittivity, whereas Wiki says material permittivity works
			{name='ionDebyeLength', value=1},
			-- normalized ion Debye length: 
			-- lambdaHat_d = lambda_d / l_r
				
			{name='min_rho', value=1e-4},
			{name='min_P', value=1e-4},
		})

		-- TODO get GridSolver's solverStruct
		-- solverVars is no longer a global variable
		-- so we have to get this from the solver
		-- but do we even need a solver_t for the base class?
		--[[
		self.solverStruct:makeType()
		self.solver_t = self.solverStruct.typename
		self.solverPtr = ffi.new(self.solver_t)
		--]]

		-- call this after we've assigned 'self' all its fields
		self:replaceSourceKernels()

		self.t = 0		
	end

	function templateClass:getConsLRTypeCode() return '' end

	function templateClass:replaceSourceKernels()
		local chargeMassRatio_ion = 1
		local chargeMassRatio_electron = .01
		local eps0 = 1 / (4 * math.pi)

		local lines = table{
			'// TwoFluidEMHDSeparateBehavior:replaceSourceKernels() end',
			self.modules:getCodeAndHeader(self.sharedModulesEnabled:keys():sort():unpack()),
--			self.eqn:getTypeCode(),

			'#define eps0 '..clnumber(eps0),
			[[
<? for _,species in ipairs{'ion', 'electron'} do 

	local speciesElectronMassRatio
	if species == 'ion' then
		speciesElectronMassRatio = 'solver->ionElectronMassRatio'
	elseif species == 'electron' then
		speciesElectronMassRatio = '1'
	end
?>

kernel void addSource_<?=species?>(
	constant <?=solver_t?> const * const solver,
	global <?=euler_cons_t?> * const derivBuf,
	global <?=euler_cons_t?> const * const UBuf,
	global <?=maxwell_cons_t?> const * const maxwellUBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	global <?=euler_cons_t?> * const deriv = derivBuf + index;
	global <?=euler_cons_t?> const * const U = UBuf + index;
	global <?=maxwell_cons_t?> const * const maxwellU = maxwellUBuf + index;
	
	deriv->m.x += (<?=speciesElectronMassRatio?> / normalizedIonLarmorRadius) * (U->rho * maxwellU->E.x + U->m.y * maxwellU->B.z - U->m.z * maxwellU->B.y);
	deriv->m.y += (<?=speciesElectronMassRatio?> / normalizedIonLarmorRadius) * (U->rho * maxwellU->E.y + U->m.z * maxwellU->B.x - U->m.x * maxwellU->B.z);
	deriv->m.z += (<?=speciesElectronMassRatio?> / normalizedIonLarmorRadius) * (U->rho * maxwellU->E.z + U->m.x * maxwellU->B.y - U->m.y * maxwellU->B.x);
	deriv->ETotal += (<?=speciesElectronMassRatio?> / normalizedIonLarmorRadius) * real3_dot(maxwellU->E, U->m);
}

<? end ?>

kernel void addSource_maxwell(
	constant <?=solver_t?> const * const solver,
	global <?=maxwell_cons_t?> * const derivBuf,
	global <?=euler_cons_t?> const * const ionUBuf,
	global <?=euler_cons_t?> const * const electronUBuf
) {
	<?=SETBOUNDS_NOGHOST?>();
	global <?=maxwell_cons_t?> * const deriv = derivBuf + index;
	global <?=euler_cons_t?> const * const ionU = ionUBuf + index;
	global <?=euler_cons_t?> const * const electronU = electronUBuf + index;
	deriv->D = real3_sub(deriv->D, real3_real_mul(ionU->m, chargeMassRatio_ion));
	deriv->D = real3_sub(deriv->D, real3_real_mul(electronU->m, chargeMassRatio_electron));
	
	
	real normalizedIonDebyeLengthSq	= normalizedIonDebyeLength * normalizedIonDebyeLength;
	deriv->E.x -= (ionU->m.x * solver->ionChargeMassRatio + elecU->m.x * elecChargeMassRatio) / (normalizedIonDebyeLengthSq * normalizedIonLarmorRadius);
	deriv->E.y -= (ionU->m.y * solver->ionChargeMassRatio + elecU->m.y * elecChargeMassRatio) / (normalizedIonDebyeLengthSq * normalizedIonLarmorRadius);
	deriv->E.z -= (ionU->m.z * solver->ionChargeMassRatio + elecU->m.z * elecChargeMassRatio) / (normalizedIonDebyeLengthSq * normalizedIonLarmorRadius);
	deriv->phi += (ionU->rho * solver->ionChargeMassRatio + elecU->rho * elecChargeMassRatio) / (normalizedIonDebyeLengthSq * normalizedIonLarmorRadius) * solver->divPhiWavespeed;
}

]],
			'// TwoFluidEMHDSeparateBehavior:replaceSourceKernels() end',
		}
		local code = lines:concat'\n'
		code = template(code, {
			solver = self,
			euler_cons_t = self.electron.eqn.symbols.cons_t,
			maxwell_cons_t = self.maxwell.eqn.symbols.cons_t,
		})
		
		-- the solver.Program is only unique to solvers wrt env and domain
		-- and those will match between these multi solvers
		-- so it doesn't matter which solver.Program I use
		self.addSourceProgramObj = self.ion.Program{name='addSource', code=code}
		self.addSourceProgramObj:compile()

		self.ion.addSourceKernelObj = self.addSourceProgramObj:kernel'addSource_ion'
		self.ion.addSourceKernelObj.obj:setArg(0, self.ion.solverBuf)
		self.ion.addSourceKernelObj.obj:setArg(2, self.ion.UBuf)
		self.ion.addSourceKernelObj.obj:setArg(3, self.maxwell.UBuf)
		
		self.electron.addSourceKernelObj = self.addSourceProgramObj:kernel'addSource_electron'
		self.electron.addSourceKernelObj.obj:setArg(0, self.electron.solverBuf)
		self.electron.addSourceKernelObj.obj:setArg(2, self.electron.UBuf)
		self.electron.addSourceKernelObj.obj:setArg(3, self.maxwell.UBuf)

		self.maxwell.addSourceKernelObj = self.addSourceProgramObj:kernel'addSource_maxwell'
		self.maxwell.addSourceKernelObj.obj:setArg(0, self.maxwell.solverBuf)
		self.maxwell.addSourceKernelObj.obj:setArg(2, self.ion.UBuf)
		self.maxwell.addSourceKernelObj.obj:setArg(3, self.electron.UBuf)
	end

	function templateClass:callAll(name, ...)
		local args = setmetatable({...}, table)
		args.n = select('#',...)
		return self.solvers:mapi(function(solver)
			return solver[name](solver, args:unpack(1, args.n))
		end):unpack()
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
			ig.igPushID_Str('subsolver '..i)
			if ig.igCollapsingHeader('sub-solver '..solver.name..':') then
				solver:updateGUI()
			end
			ig.igPopID()
		end
	end

	return templateClass
end

return TwoFluidEMHDBehavior

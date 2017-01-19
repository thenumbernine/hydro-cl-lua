--[[
using 2014 Abgrall, Kumar "Robust Finite Volume Schemes for Two-Fluid Plasma Equations

mass_ion / mass_electron is 1836 in real life
paper uses 25 
also puts ion pressure at 1/100'th the electron pressure

--]]

local class = require 'ext.class'
local table = require 'ext.table'
local ig = require 'ffi.imgui'
local vec3sz = require 'solver.vec3sz'


local TwoFluidEMHDRoe = class()

TwoFluidEMHDRoe.name = 'two-fluid EMHD Roe'

function TwoFluidEMHDRoe:init(args)
	self.app = assert(args.app)

	-- how to specify initial conditions for the both of them?
	-- separate args? maxwellInitConds vs eulerInitConds?
	-- same name in init/euler and init/maxwell?
	-- both?

	local TwoFluidIonEulerRoe = class(require 'solver.euler-roe')
	TwoFluidIonEulerRoe.prim_t = 'euler_prim_t'
	TwoFluidIonEulerRoe.cons_t = 'euler_cons_t'
	TwoFluidIonEulerRoe.consLR_t = 'euler_consLR_t'
	TwoFluidIonEulerRoe.eigen_T = 'euler_eigen_t'
	function TwoFluidIonEulerRoe:init(args)
		TwoFluidIonEulerRoe.super.init(self, args)
		self.name = 'ion '..self.name
	end
	self.ion = TwoFluidIonEulerRoe(args)

	local TwoFluidElectronEulerRoe = class(require 'solver.euler-roe')
	TwoFluidElectronEulerRoe.prim_t = 'euler_prim_t'
	TwoFluidElectronEulerRoe.cons_t = 'euler_cons_t'
	TwoFluidElectronEulerRoe.consLR_t = 'euler_consLR_t'
	TwoFluidElectronEulerRoe.eigen_T = 'euler_eigen_t'
	function TwoFluidElectronEulerRoe:init(args)
		TwoFluidElectronEulerRoe.super.init(self, args)
		self.name = 'electron '..self.name
	end
	self.electron = TwoFluidElectronEulerRoe(args)

	local TwoFluidMaxwellRoe = class(require 'solver.maxwell-roe')
	TwoFluidMaxwellRoe.prim_t = 'maxwell_prim_t'
	TwoFluidMaxwellRoe.cons_t = 'maxwell_cons_t'
	TwoFluidMaxwellRoe.consLR_t = 'maxwell_consLR_t'
	TwoFluidMaxwellRoe.eigen_T = 'maxwell_eigen_t'
	self.maxwell = TwoFluidMaxwellRoe(args)

	self.solvers = table{self.ion, self.electron, self.maxwell}

	self.displayVars = table():append(self.solvers:map(function(solver) return solver.displayVars end):unpack())

	select(2, self.ion.displayVars:find(nil, function(var) return var.name == 'U_rho' end)).enabled[0] = true 
	select(2, self.electron.displayVars:find(nil, function(var) return var.name == 'U_rho' end)).enabled[0] = true 
	select(2, self.maxwell.displayVars:find(nil, function(var) return var.name == 'U_Ex' end)).enabled[0] = false 
--	local var = select(2, self.maxwell.displayVars:find(nil, function(var) return var.name == 'U_div_B' end)) var.enabled[0] = true var.heatMapFixedRangePtr[0] = false
--	local var = select(2, self.maxwell.displayVars:find(nil, function(var) return var.name == 'U_div_E' end)) var.enabled[0] = true var.heatMapFixedRangePtr[0] = false
--	local var = select(2, self.maxwell.displayVars:find(nil, function(var) return var.name == 'ePot_0' end)) var.enabled[0] = true var.heatMapFixedRangePtr[0] = false

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
function TwoFluidEMHDRoe:getCoordMapCode() 
	return self.ion:getCoordMapCode() 
end

function TwoFluidEMHDRoe:getConsLRTypeCode() return '' end

local clnumber = require 'clnumber'
function TwoFluidEMHDRoe:replaceSourceKernels()
	local chargeMassRatio_ion = 1
	local chargeMassRatio_electron = .01
	local eps0 = 1

	require 'solver.roe'.createCodePrefix(self)

	local lines = table{
		self.codePrefix,
		self.ion.eqn:getTypeCode(),
		self.maxwell.eqn:getTypeCode(),

		'#define chargeMassRatio_ion '..clnumber(chargeMassRatio_ion),
		'#define chargeMassRatio_electron '..clnumber(chargeMassRatio_electron),
		'#define eps0 '..clnumber(eps0),
		[[

<? for _,species in ipairs{'ion', 'electron'} do ?>

kernel void addSource_<?=species?>(
	global euler_cons_t* derivBuf,
	const global euler_cons_t* UBuf,
	const global maxwell_cons_t* maxwellUBuf
) {
	SETBOUNDS(2,2);
	global euler_cons_t* deriv = derivBuf + index;
	const global euler_cons_t* U = UBuf + index;
	const global maxwell_cons_t* maxwellU = maxwellUBuf + index;
	deriv->m.x += chargeMassRatio_<?=species?> * (maxwellU->epsE.x / eps0 + U->m.y * maxwellU->B.z - U->m.z * maxwellU->B.y);
	deriv->m.y += chargeMassRatio_<?=species?> * (maxwellU->epsE.y / eps0 + U->m.z * maxwellU->B.x - U->m.x * maxwellU->B.z);
	deriv->m.z += chargeMassRatio_<?=species?> * (maxwellU->epsE.z / eps0 + U->m.x * maxwellU->B.y - U->m.y * maxwellU->B.x);
	deriv->ETotal += chargeMassRatio_<?=species?> * real3_dot(maxwellU->epsE, U->m) / eps0;
}

<? end ?>

kernel void addSource_maxwell(
	global maxwell_cons_t* derivBuf,
	const global euler_cons_t* ionUBuf,
	const global euler_cons_t* electronUBuf
) {
	SETBOUNDS(2,2);
	global maxwell_cons_t* deriv = derivBuf + index;
	const global euler_cons_t* ionU = ionUBuf + index;
	const global euler_cons_t* electronU = electronUBuf + index;
	deriv->epsE = real3_sub(deriv->epsE, real3_scale(ionU->m, chargeMassRatio_ion));
	deriv->epsE = real3_sub(deriv->epsE, real3_scale(electronU->m, chargeMassRatio_electron));
}
]]
	}
	local code = lines:concat'\n'
	code = require 'template'(code, {solver=self})
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

function TwoFluidEMHDRoe:callAll(name, ...)
	local args = setmetatable({...}, table)
	args.n = select('#',...)
	return self.solvers:map(function(solver)
		return solver[name](solver, args:unpack(1, args.n))
	end):unpack()
end

function TwoFluidEMHDRoe:getCoordMapGLSLCode()
	return self.ion:getCoordMapGLSLCode()
end

function TwoFluidEMHDRoe:createEqn()
	self:callAll'createEqn'
end

function TwoFluidEMHDRoe:resetState()
	self:callAll'resetState'
	self.t = self.ion.t
end

function TwoFluidEMHDRoe:boundary()
	self:callAll'boundary'
end

function TwoFluidEMHDRoe:calcDT()
	return math.min(self:callAll'calcDT')
end

function TwoFluidEMHDRoe:step(dt)
	self:callAll('step', dt)
	self.t = self.ion.t
end

-- same as Solver.update
-- note this means sub-solvers' update() will be skipped
-- so best to put update stuff in step()
function TwoFluidEMHDRoe:update()
	self:boundary()
	local dt = self:calcDT()
	self:step(dt)
end

function TwoFluidEMHDRoe:getTex(var) 
	return self.solverForDisplayVars[var].tex
end

function TwoFluidEMHDRoe:calcDisplayVarToTex(var)
	return self.solverForDisplayVars[var]:calcDisplayVarToTex(var)
end

function TwoFluidEMHDRoe:calcDisplayVarRange(var)
	return self.solverForDisplayVars[var]:calcDisplayVarRange(var)
end

function TwoFluidEMHDRoe:updateGUI()
	for i,solver in ipairs(self.solvers) do
		ig.igPushIdStr('subsolver '..i)
		if ig.igCollapsingHeader('sub-solver '..solver.name..':') then
			solver:updateGUI()
		end
		ig.igPopId()
	end
end

return TwoFluidEMHDRoe

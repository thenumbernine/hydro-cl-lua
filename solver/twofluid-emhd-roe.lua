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

	self.ion = require 'solver.euler-roe'(args)

	self.electron = require 'solver.euler-roe'(args)

	local maxwellArgs = table(args)
	maxwellArgs.eqn = 'maxwell'
	self.maxwell = require 'solver.roe'(maxwellArgs)

	self.solvers = table{self.ion, self.electron, self.maxwell}

	self.displayVars = table():append(self.solvers:map(function(solver) return solver.displayVars end):unpack())

	self.solverForDisplayVars = table()
	for _,solver in ipairs(self.solvers) do
		for _,var in ipairs(solver.displayVars) do
			self.solverForDisplayVars[var] = solver 
		end
	end

	
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

local clnumber = require 'clnumber'
function TwoFluidEMHDRoe:replaceSourceKernels()
	local chargeMassRatio_ion = 1
	local chargeMassRatio_electron = .01
	local eps0 = 1

	require 'solver.roe'.createCodePrefix(self)

	local lines = table{
		self.codePrefix,
		
		self.ion.eqn:getTypeCode():gsub('cons_t', 'euler_cons_t'),
		
		self.maxwell.eqn:getTypeCode():gsub('cons_t', 'maxwell_cons_t'),

		'#define chargeMassRatio_ion '..clnumber(chargeMassRatio_ion),
		'#define chargeMassRatio_electron '..clnumber(chargeMassRatio_electron),
		'#define eps0 '..clnumber(eps0),
		[[

<? for _,species in ipairs{'ion', 'electron'} do ?>

__kernel void addSource_<?=species?>(
	__global euler_cons_t* derivBuf,
	const __global euler_cons_t* UBuf,
	const __global maxwell_cons_t* maxwellUBuf
) {
	SETBOUNDS(2,2);
	__global euler_cons_t* deriv = derivBuf + index;
	const __global euler_cons_t* U = UBuf + index;
	const __global maxwell_cons_t* maxwellU = maxwellUBuf + index;
	deriv->m.x += chargeMassRatio_<?=species?> * (maxwellU->epsE.x / eps0 + U->m.y * maxwellU->B.z - U->m.z * maxwellU->B.y);
	deriv->m.y += chargeMassRatio_<?=species?> * (maxwellU->epsE.y / eps0 + U->m.z * maxwellU->B.x - U->m.x * maxwellU->B.z);
	deriv->m.z += chargeMassRatio_<?=species?> * (maxwellU->epsE.z / eps0 + U->m.x * maxwellU->B.y - U->m.y * maxwellU->B.x);
	deriv->ETotal += chargeMassRatio_<?=species?> * real3_dot(maxwellU->epsE, U->m) / eps0;
}

<? end ?>

__kernel void addSource_maxwell(
	__global maxwell_cons_t* derivBuf,
	const __global euler_cons_t* ionUBuf,
	const __global euler_cons_t* electronUBuf
) {
	SETBOUNDS(2,2);
	__global maxwell_cons_t* deriv = derivBuf + index;
	const __global euler_cons_t* ionU = ionUBuf + index;
	const __global euler_cons_t* electronU = electronUBuf + index;
	deriv->epsE = real3_sub(deriv->epsE, real3_scale(ionU->m, chargeMassRatio_ion));
	deriv->epsE = real3_sub(deriv->epsE, real3_scale(electronU->m, chargeMassRatio_electron));
}
]]
	}
	local code = lines:concat'\n'
	code = require 'processcl'(code, {solver=self})
print(require 'showcode'(code))	
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

function TwoFluidEMHDRoe:update()
	self:boundary()
	local dt = self:calcDT()
	self:step(dt)
end

function TwoFluidEMHDRoe:getTex(var) 
	return self.solverForDisplayVars[var].tex
end

function TwoFluidEMHDRoe:calcDisplayVarToTex(var)
	self.solverForDisplayVars[var]:calcDisplayVarToTex(var)
end

function TwoFluidEMHDRoe:updateGUI()
	for i,solver in ipairs(self.solvers) do
		ig.igPushIdStr('subsolver '..i)
		if ig.igCollapsingHeader('sub-solver '..solver.name..':') then
			self.ion:updateGUI()
		end
		ig.igPopId()
	end
end

return TwoFluidEMHDRoe

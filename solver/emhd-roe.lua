--[[
using 2014 Abgrall, Kumar "Robust Finite Volume Schemes for Two-Fluid Plasma Equations
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local ig = require 'ffi.imgui'
local vec3sz = require 'solver.vec3sz'

local TwoFluidEMHDRoe = class()

TwoFluidEMHDRoe.name = 'two-fluid EMHD Roe'

function TwoFluidEMHDRoe:init(args)
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

	-- change the default maxwell displayed variable
	select(2, self.maxwell.displayVars:find(nil, function(var) return var.name == 'U_Ex' end)).enabled[0] = false
	select(2, self.maxwell.displayVars:find(nil, function(var) return var.name == 'U_Ez' end)).enabled[0] = true 
	self.maxwell:refreshDisplayProgram()

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

	self.t = 0
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
		if ig.igCollapsingHeader('sub-solver '..solver.name..':') then
			ig.igPushIdStr('subsolver '..i)
			self.ion:updateGUI()
			ig.igPopId()
		end
	end
end

return TwoFluidEMHDRoe

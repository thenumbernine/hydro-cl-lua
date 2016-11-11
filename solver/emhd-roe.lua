local class = require 'ext.class'
local table = require 'ext.table'
local ig = require 'ffi.imgui'
local vec3sz = require 'solver.vec3sz'

local EMHDRoe = class()

EMHDRoe.name = 'EMHD Roe'

function EMHDRoe:init(args)
	-- how to specify initial conditions for the both of them?
	-- separate args? maxwellInitConds vs eulerInitConds?
	-- same name in init/euler and init/maxwell?
	-- both?

	self.euler = require 'solver.euler-roe'(args)
	
	local maxwellArgs = table(args)
	maxwellArgs.eqn = 'maxwell'
	self.maxwell = require 'solver.roe'(maxwellArgs)
	
	self.displayVars = table():append(self.euler.displayVars, self.maxwell.displayVars)

	self.solverForDisplayVars = table()
	for _,var in ipairs(self.euler.displayVars) do
		self.solverForDisplayVars[var] = self.euler
	end
	for _,var in ipairs(self.maxwell.displayVars) do
		self.solverForDisplayVars[var] = self.maxwell
	end

	self.numGhost = self.euler.numGhost
	
	self.color = vec3(math.random(), math.random(), math.random()):normalize()
	
	self.dim = self.euler.dim
	self.gridSize = vec3sz(self.euler.gridSize)
	self.mins = vec3(self.euler.mins:unpack())
	self.maxs = vec3(self.euler.maxs:unpack())

--[[
	self.tex = self.euler.tex
	self.maxwell.tex = self.euler.tex
	self.maxwell.texCLMem = self.euler.texCLMem
	self.maxwell.calcDisplayVarToTexPtr = self.euler.calcDisplayVarToTexPtr
--]]

	self.t = 0
end

function EMHDRoe:callAll(name, ...)
	local res1 = self.euler[name](self.euler, ...)
	local res2 = self.maxwell[name](self.maxwell, ...)
	return res1, res2
end

function EMHDRoe:getCoordMapGLSLCode()
	return self.euler:getCoordMapGLSLCode()
end

function EMHDRoe:createEqn()
	self:callAll'createEqn'
end

function EMHDRoe:resetState()
	self:callAll'resetState'
	self.t = self.euler.t
end

function EMHDRoe:boundary()
	self:callAll'boundary'
end

function EMHDRoe:calcDT()
	return math.min(self:callAll'calcDT')
end

function EMHDRoe:step(dt)
	self:callAll('step', dt)
	self.t = self.euler.t
end

function EMHDRoe:update()
	self:boundary()
	local dt = self:calcDT()
	self:step(dt)
end

function EMHDRoe:getTex(var) 
	return self.solverForDisplayVars[var].tex
end

function EMHDRoe:calcDisplayVarToTex(var)
	self.solverForDisplayVars[var]:calcDisplayVarToTex(var)
end

function EMHDRoe:updateGUI()
	if ig.igCollapsingHeader'Euler Solver:' then
		ig.igPushIdStr('euler')
		self.euler:updateGUI()
		ig.igPopId()
	end
	if ig.igCollapsingHeader'Maxwell Solver:' then
		ig.igPushIdStr('maxwell')
		self.maxwell:updateGUI()
		ig.igPopId()
	end
end

return EMHDRoe

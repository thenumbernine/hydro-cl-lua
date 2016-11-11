local class = require 'ext.class'
local ig = require 'ffi.imgui'

local EMHDRoe = class()

EMHDRoe.name = 'EMHD Roe'

function EMHDRoe:init(args)
	-- how to specify initial conditions for the both of them?
	-- separate args? maxwellInitConds vs eulerInitConds?
	-- same name in init/euler and init/maxwell?
	-- both?

	self.euler = require 'solver.euler-roe'(args)
	
	local maxwellArgs = table(args)
	maxwellArgs.eqn = require 'eqn.maxwell'()
	self.maxwell = require 'solver.roe'(args)
	
--self.numGhost = self.euler.numGhost
	self.dim = self.euler.dim
	self.gridSize = vec3sz(self.euler.gridSize)
	self.mins = vec3(self.euler.mins:unpack())
	self.maxs = vec3(self.euler.maxs:unpack())

	self.displayVars = table():append(self.euler.displayVars, self.maxwell.displayVars)

	-- honestly, does the solver need this?
	-- I think the display should have this
	-- which means put it in app
	self.tex = self.euler.tex
end

function EMHDRoe:callAll(name, ...)
	local res1 = self.euler[name](self.euler, ...)
	local res2 = self.maxwell[name](self.maxwell, ...)
	return res1, res2
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

function EMHDRoe:updateGUI()
	ig.igPushIdStr('euler')
	self.euler:updateGUI()
	ig.igPopId()
	ig.igPushIdStr('maxwell')
	self.maxwell:updateGUI()
	ig.igPopId()
end

function EMHDRoe:calcDisplayVarToTex(varIndex, var)
	-- by here varIndex will poitn to the index in the master list
	-- so we have to find it in each sub-solver's list
	
	varIndex = self.euler.displayVars:find(nil, function(var) return var.name == varName end)
	if varIndex then
		self.euler:calcDisplayVarToTex(varIndex, var)
		return
	end

	varIndex = self.maxwell.displayVars:find(nil, function(var) return var.name == varName end)
	if varIndex then
		self.maxwell:calcDisplayVarToTex(varIndex, var)
		return
	end

	error("tried to display an unknown var.  shouldn't have got here because varIndex says it's in our list")
end

return EMHDRoe

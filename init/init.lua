local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local time = table.unpack(require 'time')

--[[
name = name of the initial condition
guiVars = any gui variables that the initial conditions wants to expose
overrideGuiVars = any of the equation's gui vars that the initial conditions wants to override
initState = function(self, solver) returns the OpenCL code for the initial conditions
--]]

local InitCond = class()

function InitCond:refreshInitStateProgram(solver)
	local initStateCode = table{
		solver.codePrefix,
		self.header and self:header(solver) or '',
		solver.eqn:getInitStateCode(),
	}:concat'\n'
	
	time('compiling init state program', function()
		solver.initStateProgramObj = solver.Program{
			code = initStateCode,
		}
		solver.initStateProgramObj:compile()
	end)
	solver.initStateKernelObj = solver.initStateProgramObj:kernel('initState', solver.solverPtr, solver.UBuf)

	-- here's an ugly hack ...
	-- I need a custom init state kernel for the GLM_MHD only
	-- and it shares init conditions with a lot of other solvers
	-- so ...
	if require 'eqn.glm-mhd'.is(solver.eqn) then
		solver.initDerivsKernelObj = solver.initStateProgramObj:kernel('initDerivs', solver.UBuf)
	end
end

function InitCond:initState(solver)
	return '//no code from InitCond:initState() was provided'
end

-- called when the solver resets
function InitCond:resetState(solver)
	solver.initStateKernelObj()
	if require 'eqn.glm-mhd'.is(solver.eqn) then
		solver:boundary()
		solver.initDerivsKernelObj()
	end
	solver:boundary()
end

return InitCond

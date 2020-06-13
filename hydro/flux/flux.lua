local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local clnumber = require 'cl.obj.number'

local Flux = class()

function Flux:init(solver)
	self.solver = solver
end

function Flux:getSolverCode()
	return table{
		template(file[self.solverCodeFile], {
			solver = self.solver,
			eqn = self.solver.eqn,
			clnumber = clnumber,
		}),
	}:concat'\n'
end

return Flux

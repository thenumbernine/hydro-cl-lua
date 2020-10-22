local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local clnumber = require 'cl.obj.number'

local Flux = class()

function Flux:init(args)
	self.solver = assert(args.solver)
end

function Flux:initCodeModules()
	self.solver.modules:add{
		name = 'calcFlux',
		-- TODO these may vary depending on the solverCodeFile contents
		depends = self:getModuleDepends_calcFlux(),
		code = template(file[self.solverCodeFile], {
			solver = self.solver,
			eqn = self.solver.eqn,
			clnumber = clnumber,
		}),
	}	
end

function Flux:getModuleDepends_calcFlux()
	return {
		'solver_t',
		'cons_t',
		'waves_t',
		'normal_t',
		-- used by most subclasses:
		'cons_parallelPropagate',
		'math',	
		'solver.macros',
		'eqn.waveCode',
	}
end

return Flux

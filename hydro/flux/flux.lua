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
		name = 'Flux.calcFlux',
		-- TODO these may vary depending on the solverCodeFile contents
		depends = {
			'solver.solver_t',
			'eqn.cons_t',
			'eqn.waves_t',
			'coord.normal',
			'eqn.solvercode',
			'eqn.cons_parallelPropagate',
			
			-- specific to Roe, when solver.fluxLimiter>1
			'fluxLimiter',
			'fluxFromCons',
		},
		code = template(file[self.solverCodeFile], {
			solver = self.solver,
			eqn = self.solver.eqn,
			clnumber = clnumber,
		}),
	}	
end

return Flux

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
	self.solver.modules:addFromMarkup(
		template(file[self.solverCodeFile], {
			solver = self.solver,
			eqn = self.solver.eqn,
			clnumber = clnumber,
			flux = self,
		})
	)
end

return Flux

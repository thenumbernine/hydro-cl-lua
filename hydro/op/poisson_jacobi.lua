local table = require 'ext.table'
local class = require 'ext.class'
local path = require 'ext.path'
local template = require 'template'
local Relaxation = require 'hydro.op.relaxation'

local PoissonJacobi = class(Relaxation)

PoissonJacobi.name = 'poisson_jacobi'

PoissonJacobi.solverCodeFile = 'hydro/op/poisson.cl'

local poissonJacobiCode = path'hydro/op/poisson_jacobi.clcpp':read()

--[[
args:
	codeDepends = any additional dependencies that are used in any of the code that is provided
	i.e. maybe in a subclass (nodiv) readVectorField / writeVectorField / chargeCode / etc
--]]
function PoissonJacobi:init(args)
	PoissonJacobi.super.init(self, args)

	self.codeDepends = args.codeDepends
end

function PoissonJacobi:initCodeModules()
	local solver = self.solver
	PoissonJacobi.super.initCodeModules(self)
	solver.modules:addFromMarkup{
		code = solver.eqn:template(
			table{
				poissonJacobiCode,
				-- provided by child class:
				self:getPoissonCode() or '',
			}:concat'\n',
			table(self.symbols, {
				op = self,
			})
		),
	}
	-- this is already set in super
	solver.solverModulesEnabled[self.symbols.solveJacobi] = true
end

return PoissonJacobi

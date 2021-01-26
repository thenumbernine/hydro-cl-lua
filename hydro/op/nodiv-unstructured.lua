--[[
right now nodiv is a collection of a few classes, depending on what underlying solver you use
for the jacobi version:
nodiv.lua
	-> poisson_jacobi.lua
		-> poisson.cl (shared with poisson_krylov)
		-> poisson_jacobi.cl
		-> relaxation.lua

for the krylov version:
nodiv.lua
	-> poisson.cl (shared with poisson_jacobi)
	-> solver/cl/$solverClass.lua

what do the classes fulfill externally?
- their :step(dt) function removes divergence from the field specified in cmdline args.

what goes on within the jacobi version?
NoDiv:step()
	call NoDiv:relax()
		Relaxation:relax()
		call solveJacobiKernelObj()
		call copyWriteToPotentialNoGhsotKernelObj
		call self:potentialBoundary()
		do stop on epsilon stuff maybe
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local CLBuffer = require 'cl.obj.buffer'

function NoDiv:resetState()
	local solver = self.solver
	if self.enableField and not solver[self.enableField] then return end
	self.initPotentialKernelObj()
	self:potentialBoundary()
	self:relax()
end

function NoDiv:init(args)
	local solver = assert(args.solver)
	self.solver = solver
	
	require 'hydro.code.symbols'(self, self:getSymbolFields())

	self.writeBufObj = CLBuffer{
		env = solver.app.env,
		name = 'writeBuf',
		type = solver.app.real,
		count = solver.numCells,
	}
end

function NoDiv:getSymbolFields()
	return table{
		'initPotential',
	}
end

function NoDiv:initCodeModules()
	local solver = self.solver
	solver.modules:add{
		name = name,
		depends = {
			solver.symbols.SETBOUNDS_NOGHOST,
		},
		code = solver.eqn:template(file[self.solverCodeFile], {op = self}),
	}
	solver.solverModulesEnabled[name] = true
end



return NoDiv

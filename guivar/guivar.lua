local class = require 'ext.class'

local GuiVar = class()

function GuiVar:init(args)
	self.name = assert(args.name)
	self.onChange = args.onChange
	
	-- set to true or the var will be a runtime var
	--  this kind of var inserts itself into solver_t
	-- don't change this at runtime
	self.compileTime = args.compileTime
end

function GuiVar:refresh(value, solver)
	print('refreshing '..self.name..' = '..tostring(value))
	if self.onChange then self:onChange(value, solver) end
	if self.compileTime then
		solver:refreshCodePrefix()
	else
		self:setToSolver(solver)
		solver:refreshSolverBuf()
	end
end

return GuiVar

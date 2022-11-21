local class = require 'ext.class'

local GuiVar = class()

function GuiVar:init(args)
	self.name = assert(args.name)
	self.onChange = args.onChange
	self.units = args.units

	-- set to true or the var will be a runtime var
	--  this kind of var inserts itself into solver_t
	-- don't change this at runtime
	self.compileTime = args.compileTime
end

function GuiVar:refresh(solver)
	local value = self:getValue()
	print('refreshing '..self.name..' = '..tostring(value))
	if self.onChange then self:onChange(value, solver) end
	if self.compileTime then
		solver:refreshCodePrefix()
	else
		self:refreshInStruct(solver)
	end
end

--[[
TODO think of a better name
this gets .value except for combo then it gets the .options[] of .value
--]]
function GuiVar:getValue()
	return self.value
end

-- what struct within solver to read/write our values from?
-- default is the solverPtr
-- however initCond guiVars will want to read/write initCondPtr
function GuiVar:refreshInStruct(solver)
	-- use self.value instead of 'value' because in combo 'self.value' is the key and 'value == :getValue()' is the key's value / desc of the key
	solver.solverPtr[self.name] = self.value
	solver:refreshSolverBuf()
end

return GuiVar

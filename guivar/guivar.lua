local class = require 'ext.class'

local GuiVar = class()

function GuiVar:init(args)
	self.name = assert(args.name)
end

function GuiVar:refresh(value, solver)
	print('refreshing '..self.name..' = '..value)
	
	-- this in turn calls eqn:getCodePrefix, which gets the guivar code
	solver:createCodePrefix()
	solver:refreshInitStateProgram()
	solver:refreshSolverProgram()
	solver:refreshDisplayProgram()
end

return GuiVar

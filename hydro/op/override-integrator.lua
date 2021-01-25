--[[
args:
	op = which op class to override
		should use addSource to add to the solver
	integrator = which integrator class to override it with
--]]
return function(args)
	local class = require 'ext.class'

	local opClass = assert(args.op)
	local integratorClass = assert(args.integrator)

	local behavior = class(opClass) 

	function behavior:refreshSolverProgram()
		if opClass.refreshSolverProgram then
			opClass.refreshSolverProgram(self)
		end

		self.integrator = integratorClass(self.solver)
	end

	-- disable the parent class' 'addSource'
	behavior.addSource = nil

	-- instead replace it with an 'update'
	-- and do the 'addSource' manually here

	function behavior:step(dt)
		self.integrator:integrate(dt, function(derivBufObj)
			opClass.addSource(self, derivBufObj)
		end)

		if opClass.step then
			opClass.step(self, dt)
		end
	end

	return behavior
end

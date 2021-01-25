-- TODO this whole thing could be a behavior for transforming any op's addSource explicit update into an implicit, or any update using an alternative integrator to the solver's default integrator

local class = require 'ext.class'
local ViscousExplicit = require 'hydro.op.viscous-explicit'

local ViscousImplicit = class(ViscousExplicit)

ViscousImplicit.name = 'ViscousImplicit'

function ViscousImplicit:refreshSolverProgram()
	ViscousImplicit.super.refreshSolverProgram(self)

	self.integrator = require 'hydro.int.be'(self.solver)
end

-- disable the parent class' 'addSource'
local oldAddSource = ViscousImplicit.addSource
ViscousImplicit.addSource = nil

-- instead replace it with an 'update'
-- and do the 'addSource' manually here

function ViscousImplicit:step(dt)
	self.integrator:integrate(dt, function(derivBufObj)
		oldAddSource(self, derivBufObj)
	end)
end

return ViscousImplicit 

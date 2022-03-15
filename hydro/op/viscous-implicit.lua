local ViscousImplicit = require 'hydro.op.override-integrator'{
	op = require 'hydro.op.viscous-explicit',
	integrator = require 'hydro.int.be'
}

ViscousImplicit.name = 'ViscousImplicit'

return ViscousImplicit 

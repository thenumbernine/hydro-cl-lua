--[[
me trying to do twofluid-emhd as a composite of two euler's and glm-maxwell
	ion = euler
	elec = euler
	em = glm-maxwell

what do I have to change to get it to match the manually create twofluid-emhd.lua

* would be nice to manually name the eqn1 eqn2 eqn3 fields to things like 'ion', 'elec', 'em'
* the ePot fields from both ion and elec need to be merged.
* sigma field doesn't exist.  I think it's a global constant, probably within Larmor gyroradius maybe?
* sqrt_1_eps and sqrt_1_mu don't exist.  They are assumed to be vacuum constants.
* for init cond, each euler's vars are assigned with scalars of the init cond provided vars ...
	for {rho, v, P}
	ion's scalars are {1, 1, 1/ionElectronMassRatio}
	elec's scalars are {1/ionElectronMassRatio, 1, 5}
* ops (esp NoDiv and SelfGrav) need to work.
* display vars applied per sub eqn would be nice.
* addSource needs to be able to add extra terms:
	Lorentz force:
		ion_m += ionChargeMassRatio / unit_C_per_kg * (rho * D / eps + ion_m cross B)
		ion_ETotal += ionChargeMassRatio / unit_C_per_kg * D / eps dot ion_m
		elec_m -= elecChargeMassRatio / unit_C_per_kg * (elec_rho * D / eps + elec_m cross B)
		elec_ETotal -= elecChargeMassRatio / unit_C_per_kg * D / eps dot elec_m
	Maxwell source:
		D -= (ion_m * ionChargeMassRatio + elec_m * elecChargeMassRatio) / unit_C_per_kg
		phi += eps * (ion_rho * ionChargeMassRatio + elec_rho * elecChargeMassRatio) / unit_C_per_kg * divPhiWavespeed / unit_m_per_s
	...and any connection coefficients ... which should already be within the previous subeqns addSource
	... so which of these is not already in a subeqn addSource?
	right now glm-maxwell addSource has:
		D -= D / eps * sigma ... which we don't want.  this is replaced with the ion_m + elec_m above.
		and then some gradient smoothing stuff ... which we do want
		and then some connection stuff ... which we do want
	right now euler addSource has:
		connection stuff (disabled?) ... which we want
		so we would have to add the Lorentz stuff
	* addSource => addSourceCell, just like calcDT and applyInitCond

... this makes me think, for glm-maxwell, for maxwell, to have an option to either provide ...
	* eps as a field or a constant (current glm-maxwell is a field, twofluid-emhd is a constant)
	* mu as a field or a constant (current glm-maxwell is a field, twofluid-emhd is a constant)
	* sigma as a field?  at all?  what about nonlinear responses?  as a function?
--]]

local Composite = require 'hydro.eqn.composite'

local TwoFluidEMHDComposite = class(Composite)

TwoFluidEMHDComposite.name = 'twofluid_emhd_composite'

function TwoFluidEMHDComposite:init(args)
	args.subeqns = {
		{name='euler', field='ion'},
		{name='euler', field='elec'},
		{name='maxwell', field='em'},
	}
end

return TwoFluidEMHDComposite 

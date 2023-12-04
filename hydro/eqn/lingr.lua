--[[
See my "Einstein Equations, Weak Field, De-Donder Gauge" worksheet for a start of this.  or just look at gravito-electromagnetics.

This is General Relativity with weak field approximation, de-Donder gauge, and a cartesian metric (like the twofluid_plasma_lingr file),
but unlike the twofluid_plasma_lingr file, this is without the h_ab,tt = 0 approximation, therefore it is not a gravito-electromagnetism model.
Also it is only valid as cartesian vector components.

Actually same with twofluid-emhd-lingr's derivation -- only Cartesian -- but people still take it and recast the laplacian into non-Cartesian coordinates (is that still valid?)

TODO call this wave-gr instead of lin-gr?  because this is different from the GEM which assumes h_ab,tt = 0

This means when the metric perturbation tensors are inserted into the Einstein field equations they become the 10 wave equations:

h_ab,tt - h_ab,ii = 16 pi (T_ab + 1/2 g_ab (T_tt - T_ii))

and this means that the only feedback from the wave equations happens along the trace.
this then becomes ...

Pi_ab = h_ab,t
Psi_ab = h_ab,i

Pi_ab,t - Psi_ab,i = ...
Pi_ab,i + Psi_ab,i = 0

however we have some options for stress-energies

let's do ordinary euler fluid equations first

TODO how about a composite 'equation' ? instead of a composite solver file?
--]]
local Equation = require 'hydro.eqn.eqn'

local LinGR = Equation:subclass()
LinGR.name = 'lingr'
LinGR.usePressure = false

LinGR.initConds = require 'hydro.init.euler':getList()	 -- use rho as our initial condition

function LinGR:init(args)
	-- all unitless, all lower
	-- TODO just make a sym4 struct?
	self.consVars = {
		{name='Pi_tt', type='real', variance=''},		-- Pi_tt = h_tt,t
		{name='Pi_ti', type='real3', variance='l'},		-- Pi_ti.s[i] = h_ti,t
		{name='Pi_ij', type='real3s3', variance='ll'},		-- Pi_ij.s[sym(i,j)] = h_ij,t
		
		-- derivative in the prefix
		{name='Psi_ttk', type='real3', variance='l'},		-- Psi_ttk.s[k] = h_tt,k
		{name='Psi_tik', type='real3x3', variance='ll'},	-- Psi_tik.s[k].s[i] = h_ti,k
		{name='Psi_ijk', type='real3x3s3', variance='lll'},	-- Psi_ijk.s[k].s[sym(i,j)] = h_ij,k
	}

	LinGR.super.init(self, args)

	if not (
		self.solver.coord.vectorComponent == 'cartesian'
		or self.solver.coord.name == 'cartesian'
	) then
		print("This solver was derived specifically for cartesian coord.vectorComponent.  Using another coord.vectorComponent would introduce a lot more terms (TODO do this).")
		print("You are using coord.vectorComponent: "..self.solver.coord.vectorComponent)
	end
end

function LinGR:createInitState()
	LinGR.super.createInitState(self)
	self:addGuiVars{
		{name='wavespeed', value=1, units='m/s'},
	}
end


LinGR.solverCodeFile = 'hydro/eqn/lingr.cl'

-- don't use default
function LinGR:initCodeModule_fluxFromCons() end

function LinGR:eigenWaveCodePrefix(args)
	return self:template([[
real wavespeed = solver->wavespeed / unit_m_per_s;
real nLen = normal_len(<?=n?>);
]], args)
end

-- safe to do as long as eigenWaveCodePrefix doesn't use "eig"
LinGR.consWaveCodePrefix = LinGR.eigenWaveCodePrefix

function LinGR:consWaveCode(args)
	local waveIndex = args.waveIndex
	waveIndex = waveIndex % 4
	if waveIndex == 0 then
		return '-wavespeed * nLen'
	elseif waveIndex == 1 or waveIndex == 2 then
		return '0.'
	elseif waveIndex == 3 then
		return 'wavespeed * nLen'
	end
	error'got a bad waveIndex'
end
LinGR.eigenWaveCode = LinGR.consWaveCode

return LinGR

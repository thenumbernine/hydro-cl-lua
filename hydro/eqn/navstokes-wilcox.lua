--[[
Taken from Petrova - Finite Volume Methods: Powerful Means of Engineering Design
cited as 1998 Wilcox' version of Navier-Stokes

k-omega turbulence model of Navier Stokes method for finite volume
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local materials = require 'hydro.materials'
local Equation = require 'hydro.eqn.eqn'

local NavierStokesWilcox = class(Equation)
NavierStokesWilcox.name = 'navstokes_wilcox'

NavierStokesWilcox.numWaves = 7	-- v-a, v,v,v,v,v, v+a
NavierStokesWilcox.numIntStates = 7

NavierStokesWilcox.roeUseFluxFromCons = true

NavierStokesWilcox.initConds = require 'hydro.init.euler':getList()

function NavierStokesWilcox:init(args)
	self.primVars = table{
		{name='rhoBar', type='real'},
		{name='vTilde', type='real3', variance='u'},
		{name='PStar', type='real'},
		{name='k', type='real'},
		{name='omega', type='real'},
		{name='ePot', type='real'},
	}

	self.consVars = table{
		{name='rhoBar', type='real'},
		{name='rhoBar_vTilde', type='real3', variance='u'},
		{name='rhoBar_eTotalTilde', type='real', variance=''},	-- manually specify variance to variables that use underscores
		{name='rhoBar_k', type='real', variance=''},
		{name='rhoBar_omega', type='real', variance=''},
		{name='ePot', type='real'},
	}

	if args.incompressible then
		self.consVars:insert{name='mPot', type='real', units='kg/(m*s)'}
		self.primVars:insert{name='mPot', type='real', units='kg/(m*s)'}
	end

	NavierStokesWilcox.super.init(self, args)


	if require 'hydro.solver.meshsolver'.is(self.solver) then
		print("not using selfgrav with mesh solvers yet")
	else
--[[ might have to change the selfgrav terms ...
		local SelfGrav = require 'hydro.op.selfgrav'
		self.gravOp = SelfGrav{solver = self.solver}
		self.solver.ops:insert(self.gravOp)
--]]
		if args.incompressible then
			local NoDiv = require 'hydro.op.nodiv'{
				poissonSolver = require 'hydro.op.poisson_jacobi',	-- krylov is having errors.  TODO bug in its boundary code?
			}
			self.solver.ops:insert(NoDiv{
				solver = self.solver,
				vectorField = 'rhoBar_vTilde',
				potentialField = 'mPot',
			
				-- div v = 0
				-- div (m/ρ) = 0
				-- 1/ρ div m - 1/ρ^2 m dot grad ρ = 0
				-- div m = (m dot grad ρ)/ρ 
				chargeCode = self:template[[
	<? for j=0,solver.dim-1 do ?>{
		global <?=cons_t?> const * const Ujm = U - solver->stepsize.s<?=j?>;
		global <?=cons_t?> const * const Ujp = U + solver->stepsize.s<?=j?>;
		real drho_dx = (Ujp->rhoBar - Ujm->rhoBar) * (.5 / solver->grid_dx.s<?=j?>);
		source -= drho_dx * U->rhoBar_vTilde.s<?=j?> / U->rhoBar;
	}<? end ?>
]],
			})
		end
	end
end

function NavierStokesWilcox:createInitState()
	NavierStokesWilcox.super.createInitState(self)
--in order to make things work, gamma needs to be set *HERE AND IN INIT/EULER*
-- which means it is being read and written in multiple places
-- TODO consolidate that
	self:addGuiVars{
		{name='C_v', value=materials.Air.C_v},
		{name='C_p', value=materials.Air.C_p},

		-- specific gas constant R_spec
		{name='gasConstant', value=(materials.Air.C_p - materials.Air.C_v)},
		
		{name='heatCapacityRatio', value=materials.Air.C_p / materials.Air.C_v},
	}
end

-- don't use default
function NavierStokesWilcox:initCodeModule_consFromPrim_primFromCons() end
function NavierStokesWilcox:initCodeModule_fluxFromCons() end

function NavierStokesWilcox:getModuleDepends_waveCode() 
	return {
		'normal_t',
		self.symbols.eqn_common,
	}
end

function NavierStokesWilcox:getModuleDepends_displayCode() 
	return {
		self.symbols.eqn_common,
	}
end

NavierStokesWilcox.solverCodeFile = 'hydro/eqn/navstokes-wilcox.cl'

NavierStokesWilcox.displayVarCodeUsesPrims = true

function NavierStokesWilcox:getDisplayVars()
	local vars = NavierStokesWilcox.super.getDisplayVars(self)
	vars:append{
		{name='vTilde', code='value.vreal3 = W.vTilde;', type='real3'},
		{name='PStar', code='value.vreal = W.PStar;'},
		{name='eIntTilde', code='value.vreal = calc_eIntTilde(solver, &W);'},
		{name='eKinTilde', code='value.vreal = calc_eKinTilde(&W, x);'},
		{name='eTotalTilde', code='value.vreal = U->rhoBar_eTotalTilde / W.rhoBar;'},
		{name='EIntTilde', code='value.vreal = calc_EIntTilde(solver, &W);'},
		{name='EKinTilde', code='value.vreal = calc_EKinTilde(&W, x);'},
		{name='EPot', code='value.vreal = U->rhoBar * U->ePot;'},
		{name='S', code='value.vreal = W.PStar / pow(W.rhoBar, R_over_C_v + 1. );'},
		--{name='H', code='value.vreal = calc_H(W.PStar);'},
		--{name='h', code='value.vreal = calc_h(W.rhoBar, W.PStar);'},
		--{name='HTotal', code='value.vreal = calc_HTotal(W.PStar, U->rhoBar_eTotalTilde);'},
		--{name='hTotal', code='value.vreal = calc_hTotal(W.rhoBar, W.PStar, U->rhoBar_eTotalTilde);'},
		{name='Speed of Sound', code='value.vreal = calc_Cs(solver, &W);'},
		--{name='Mach number', code='value.vreal = coordLen(W.vTilde, x) / calc_Cs(solver, &W);'},
	}:append{self.gravOp and
		{name='gravity', code=self:template[[
	if (OOB(1,1)) {
		value.vreal = 0.;
	} else {
		<? 
for side=0,solver.dim-1 do ?>{
			global const <?=cons_t?>* Um = U - solver->stepsize.s<?=side?>;
			global const <?=cons_t?>* Up = U + solver->stepsize.s<?=side?>;
			value_real3->s<?=side?> = -(Up-><?=eqn.gravOp.potentialField?> - Um-><?=eqn.gravOp.potentialField?>) / (2. * cell_dx<?=side?>(x));
		}<? 
end
for side=solver.dim,2 do ?>
		value_real3->s<?=side?> = 0.;
<? end ?>
	}
]], type='real3'} or nil
	}:append{
		{name='temp', code='value.vreal = calc_eIntTilde(solver, &W) / solver->C_v;'},
	}

	vars:insert(self:createDivDisplayVar{
		field = 'vTilde', 
		getField = function(U, j)
			return U..'->rhoBar_vTilde.s'..j..' / '..U..'->rhoBar'
		end,
		units = 'kg/(m^3*s)',
	})

	vars:insert(self:createCurlDisplayVar{
		field = 'vTilde',
		getField = function(U, j)
			return U..'->rhoBar_vTilde.s'..j..' / '..U..'->rhoBar'
		end,
		units = 'm/s^2',
	})


	return vars
end

NavierStokesWilcox.eigenVars = table{
	-- Roe-averaged vars
	{name='rhoBar', type='real'},
	{name='vTilde', type='real3'},
	{name='hTotal', type='real'},
	{name='k', type='real'},
	{name='omega', type='real'},
	-- derived vars
	{name='vTildeSq', type='real'},
	{name='Cs', type='real'},
}

function NavierStokesWilcox:eigenWaveCodePrefix(n, eig, x)
	return self:template([[
real Cs_nLen = <?=eig?>->Cs * normal_len(<?=n?>);
real v_n = normal_vecDotN1(<?=n?>, <?=eig?>->vTilde);
]], {
		n = n,
		eig = '('..eig..')',
		x = x,
	})
end

function NavierStokesWilcox:consWaveCodePrefix(n, U, x)
	return self:template([[
prim_t W;
<?=primFromCons?>(&W, solver, <?=U?>, <?=x?>);
real Cs_nLen = calc_Cs(solver, &W) * normal_len(<?=n?>);
real v_n = normal_vecDotN1(<?=n?>, W.vTilde);
]], {
		n = n,
		U = '('..U..')',
		x = x,
	})
end

function NavierStokesWilcox:consWaveCode(n, eig, x, waveIndex)
	if waveIndex == 0 then
		return '(v_n - Cs_nLen)'
	elseif waveIndex >= 1 and waveIndex <= 5 then
		return 'v_n'
	elseif waveIndex == 6 then
		return '(v_n + Cs_nLen)'
	end
	error'got a bad waveIndex'
end

NavierStokesWilcox.eigenWaveCode = NavierStokesWilcox.consWaveCode

return NavierStokesWilcox

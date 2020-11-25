--[[
Stone et al 2008 - https://arxiv.org/pdf/0804.0402v1.pdf
based on Athena's version of eigenvectors of derivative of adiabatic MHD flux wrt primitives
ideal-mhd, divergence-free, conservative-based eigensystem
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local constants = require 'hydro.constants'
local Struct = require 'hydro.code.struct'
local Equation = require 'hydro.eqn.eqn'

local MHD = class(Equation)

MHD.name = 'MHD'

MHD.numWaves = 7
MHD.numIntStates = 8

MHD.roeUseFluxFromCons = true

-- hmm, we want init.euler and init.mhd here ...
MHD.initConds = require 'hydro.init.euler':getList()


-- these are calculated based on cell-centered (or extrapolated) conserved vars
-- they are used to calculate the eigensystem at a cell center or edge 
MHD.roeVars = table{
	{name='rho', type='real'},
	{name='v', type='real3'},
	{name='hTotal', type='real'},
	{name='B', type='real3'},
	{name='X', type='real'},
	{name='Y', type='real'},
}

-- here's the variables that an eigensystem uses to compute a left, right, or flux transform 
MHD.eigenVars = table(MHD.roeVars):append{

	{name='hHydro', type='real'},
	{name='aTildeSq', type='real'},

	{name='Cs', type='real'},
	{name='CAx', type='real'},
	{name='Cf', type='real'},

	{name='BStarPerpLen', type='real'},
	{name='betaY', type='real'},
	{name='betaZ', type='real'},
	{name='betaStarY', type='real'},
	{name='betaStarZ', type='real'},
	{name='betaStarSq', type='real'},

	{name='alphaF', type='real'},
	{name='alphaS', type='real'},

	{name='sqrtRho', type='real'},
	{name='sbx', type='real'},
	{name='Qf', type='real'},
	{name='Qs', type='real'},
	{name='Af', type='real'},
	{name='As', type='real'},
}

function MHD:init(args)
	
	-- TODO redo the mhd equations for a background grid metric, and take note of covariance/contravariance
	self.primVars = table{
		{name='rho', type='real', units='kg/m^3'},
		{name='v', type='real3', units='m/s', variance='u'},
		{name='P', type='real', units='kg/(m*s^2)'},
		{name='B', type='real3', units='kg/(C*s)', variance='l'},
		{name='psi', type='real', units='kg/(C*s)'},
		{name='ePot', type='real', units='m^2/s^2'},
	}

	self.consVars = table{
		{name='rho', type='real', units='kg/m^3'},
		{name='m', type='real3', units='kg/(m^2*s)', variance='u'},
		{name='ETotal', type='real', units='kg/(m*s^2)'},
		{name='B', type='real3', units='kg/(C*s)', variance='l'},
		{name='psi', type='real', units='kg/(C*s)'},
		{name='ePot', type='real', units='m^2/s^2'},
	}
	
	if args.incompressible then
		self.consVars:insert{name='mPot', type='real', units='kg/(m*s)'}
		self.primVars:insert{name='mPot', type='real', units='kg/(m*s)'}
	end

	MHD.super.init(self, args)
	
	local solver = self.solver

	self.roeStruct = Struct{solver=solver, name='roe_t', vars=self.roeVars}
	self.roeStruct:makeType()
	self.roe_t = self.roeStruct.typename

	if require 'hydro.solver.meshsolver'.is(solver) then
		print("not using ops (selfgrav, nodiv, etc) with mesh solvers yet")
	else
		if solver.dim > 1 then
			local NoDiv = require 'hydro.op.nodiv'()
			solver.ops:insert(NoDiv{
				solver = solver,
				potentialField = 'psi',
			})
		end

		local SelfGrav = require 'hydro.op.selfgrav'
		self.gravOp = SelfGrav{solver=solver}
		solver.ops:insert(self.gravOp)
	
		if args.incompressible then
			local NoDiv = require 'hydro.op.nodiv'{
				poissonSolver = require 'hydro.op.poisson_jacobi',	-- krylov is having errors.  TODO bug in its boundary code?
			}
			self.solver.ops:insert(NoDiv{
				solver = self.solver,
				vectorField = 'm',
				potentialField = 'mPot',
			
				-- div v = 0
				-- div (m/ρ) = 0
				-- 1/ρ div m - 1/ρ^2 m dot grad ρ = 0
				-- div m = (m dot grad ρ)/ρ 
				chargeCode = self:template[[
	<? for j=0,solver.dim-1 do ?>{
		global <?=eqn.cons_t?> const * const Ujm = U - solver->stepsize.s<?=j?>;
		global <?=eqn.cons_t?> const * const Ujp = U + solver->stepsize.s<?=j?>;
		real const drho_dx = (Ujp->rho - Ujm->rho) * (.5 / solver->grid_dx.s<?=j?>);
		source -= drho_dx * U->m.s<?=j?> / U->rho;
	}<? end ?>
]],
			})
		end
	end
end

--[[
B^2 = kg^2/(C^2 s^2)
B^2 / mu = kg^2/(C^2 s^2) * C^2/(kg*m) = kg/(m s^2)
--]]
MHD.guiVars = {
	{name='heatCapacityRatio', value=2},	-- 5/3 for most problems, but 2 for Brio-Wu, so I will just set it here for now (in case something else is reading it before it is set there)
	
	-- works good as mu0 = 1
	-- solver->mu0 / unit_kg_m_per_C2 == 1
	-- solver->mu0 / (unit_kg * unit_m / (unit_C * unit_C)) == 1
	-- unit_C = sqrt((unit_kg * unit_m) / solver->mu0);
	{name='mu0', value=constants.vacuumPermeability_in_kg_m_per_C2, units='(kg*m)/C^2'},
}

function MHD:initCodeModules()
	MHD.super.initCodeModules(self)
	
	self.solver.modules:add{
		name = 'roe_t',
		structs = {self.roeStruct},
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.roe_t..' roe_t;',
	}

	-- where do I put this to make it the default value for MHD solvers,
	-- but not override a value set by the init state?
	-- this used to be in Equation:getInitCondCode
	self.guiVars.coulomb.value = math.sqrt(self.guiVars.kilogram.value * self.guiVars.meter.value / self.guiVars.mu0.value)
end

-- don't use default
function MHD:initCodeModule_fluxFromCons() end
function MHD:initCodeModulePrimCons() end

function MHD:getModuleDepends_waveCode() 
	return {
		'calcCellMinMaxEigenvalues',
	}
end

function MHD:getModuleDepends_displayCode() 
	return {
		'eqn.common',
		'consFromPrim',
	}
end

MHD.solverCodeFile = 'hydro/eqn/mhd.cl'

MHD.displayVarCodeUsesPrims = true

MHD.predefinedDisplayVars = {
	'U rho',
	'U m',
	'U ETotal',
	'U D',
	'U B',
	'U div B',
}

function MHD:getDisplayVars()
	local vars = MHD.super.getDisplayVars(self)
	vars:append{
		{name='v', code='value.vreal3 = W.v;', type='real3', units='m/s'},
		{name='P', code='value.vreal = W.P;', units='kg/(m*s^2)'},
		{name='PMag', code='value.vreal = calc_PMag(solver, &W, x);', units='kg/(m*s^2)'},
		{name='PTotal', code='value.vreal = W.P + calc_PMag(solver, &W, x);', units='kg/(m*s^2)'},
		{name='eInt', code='value.vreal = calc_eInt(solver, &W);', units='m^2/s^2'},
		{name='EInt', code='value.vreal = calc_EInt(solver, &W);', units='kg/(m*s^2)'},
		{name='eKin', code='value.vreal = calc_eKin(&W, x);', units='m^2/s^2'},
		{name='EKin', code='value.vreal = calc_EKin(&W, x);', units='kg/(m*s^2)'},
		{name='eHydro', code='value.vreal = calc_eHydro(solver, &W, x);', units='m^2/s^2'},
		{name='EHydro', code='value.vreal = calc_EHydro(solver, &W, x);', units='kg/(m*s^2)'},
		{name='EM energy', code='value.vreal = calc_EM_energy(solver, &W, x);', units='kg/(m*s^2)'},
		{name='eTotal', code='value.vreal = U->ETotal / W.rho;', units='m^2/s^2'},
		{name='S', code='value.vreal = W.P / pow(W.rho, (real)solver->heatCapacityRatio);'},
		{name='H', code='value.vreal = calc_H(solver, W.P);', units='kg/(m*s^2)'},
		{name='h', code='value.vreal = calc_H(solver, W.P) / W.rho;', units='m^2/s^2'},
		{name='HTotal', code='value.vreal = calc_HTotal(solver, &W, U->ETotal, x);', units='kg/(m*s^2)'},
		{name='hTotal', code='value.vreal = calc_hTotal(solver, &W, U->ETotal, x);', units='m^2/s^2'},
		{name='speed of sound', code='value.vreal = calc_Cs(solver, &W);', units='m/s'},
		{name='alfven velocity', code='value.vreal3 = calc_CA(solver, &U);', type='real3', units='m/s'},
		{name='Mach number', code='value.vreal = coordLen(W.v, x) / calc_Cs(solver, &W);'},
		{name='temperature', code=self:template[[
<? local clnumber = require 'cl.obj.number' ?>
<? local materials = require 'hydro.materials' ?>
#define C_v				<?=('%.50f'):format(materials.Air.C_v)?>
	value.vreal = calc_eInt(solver, &W) / C_v;
]], units='K'},
		{name='primitive reconstruction error', code=self:template[[
	//prim have just been reconstructed from cons
	//so reconstruct cons from prims again and calculate the difference
	<?=eqn.cons_t?> U2;
	consFromPrim(&U2, solver, &W, x);
	value.vreal = 0;
	for (int j = 0; j < numIntStates; ++j) {
		value.vreal += fabs(U->ptr[j] - U2.ptr[j]);
	}
]]},
	}

	if self.gravOp then
		vars:insert{
			name='gravity', 
			code=self:template[[
	if (!OOB(1,1)) {
		value.vreal3 = calcGravityAccel<?=eqn.gravOp.name?>(solver, U);
	}
]], 
			type='real3', 
			units='m/s^2',
		}
	end

	vars:insert(self:createDivDisplayVar{
		field = 'v', 
		getField = function(U, j)
			return U..'->m.s'..j..' / '..U..'->rho'
		end,
		units = '1/s',
	} or nil)

	vars:insert(self:createCurlDisplayVar{
		field = 'v',
		getField = function(U, j)
			return U..'->m.s'..j..' / '..U..'->rho'
		end,
		units = '1/s',
	} or nil)

	vars:insert(self:createDivDisplayVar{field='B', units='kg/(C*m*s)'} or nil)
	vars:insert(self:createCurlDisplayVar{field='B', units='kg/(C*m*s)'} or nil)

	return vars
end

function MHD:eigenWaveCode(n, eig, x, waveIndex)
	eig = '('..eig..')'
	return ({
		eig..'->v.x - '..eig..'->Cf',
		eig..'->v.x - '..eig..'->CAx',
		eig..'->v.x - '..eig..'->Cs',
		eig..'->v.x',
		eig..'->v.x + '..eig..'->Cs',
		eig..'->v.x + '..eig..'->CAx',
		eig..'->v.x + '..eig..'->Cf',
	})[waveIndex+1] or error("got a bad waveIndex")
end

function MHD:consWaveCodePrefix(n, U, x)
	return self:template([[
range_t lambda;
calcCellMinMaxEigenvalues(&lambda, solver, <?=U?>, <?=x?>, <?=n?>); 
]], {
		n = n,
		U = '('..U..')',
		x = x,
	})
end

function MHD:consMinWaveCode(side, U, x)
	return 'lambda.min'
end

function MHD:consMaxWaveCode(side, U, x)
	return 'lambda.max'
end

return MHD

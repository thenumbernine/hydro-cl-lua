--[[
started with the mhd.lua and mhd.cl
tweaked it while looking at
2006 Xueshang et al - A 3rd Order WENO GLM-MHD Scheme for Magnetic Reconnection
2009 Mignone, Tzeferacos - A Second-Order Unsplit Godunov Scheme for Cell-Centered MHD- the CTU-GLM scheme
--]]

local table = require 'ext.table'
local constants = require 'hydro.constants'
local Struct = require 'struct'
local Equation = require 'hydro.eqn.eqn'

local MHD = Equation:subclass()

MHD.name = 'glm_mhd'

MHD.numWaves = 9
MHD.numIntStates = 9

MHD.roeUseFluxFromCons = true
MHD.useConstrainU = true

-- TODO this is broken
MHD.useFixedCh = true -- true = use a gui var, false = calculate by max(|v_i|+Cf)


MHD.use2009MignoneEqn3 = false	-- 2003 Mignone eqn 3
MHD.use2002DednerEqn24 = true	-- 2002 Dedner eqn 24, EGLM
MHD.use2002DednerEqn38 = true	-- 2002 Dedner eqn 38

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
	{name='Ch', type='real'},
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

	local solver = assert(args.solver)
	if require 'hydro.solver.meshsolver':isa(solver) then
		print("not divergence with mesh solvers yet")
		-- these don't work with mesh solver
		self.use2002DednerEqn24 = false
		self.use2002DednerEqn38 = false
	end

	MHD.super.init(self, args)
	
	self.symbols.roe_t = solver.app:uniqueName'roe_t'
	self.roeStruct = Struct{
		name = self.symbols.roe_t,
		cdef = false,
		fields = self.roeVars,
	}.class

	local UpdatePsi = require 'hydro.op.glm-mhd-update-psi'
	solver.ops:insert(UpdatePsi{solver=solver})
	
	if require 'hydro.solver.meshsolver':isa(solver) then
		print("not using ops (selfgrav, nodiv, etc) with mesh solvers yet")
	else
		local SelfGrav = require 'hydro.op.selfgrav'
		self.gravOp = SelfGrav{solver=solver}
		solver.ops:insert(self.gravOp)
	
		if args.incompressible then
			local NoDiv = require 'hydro.op.nodiv'{
				poissonSolver = require 'hydro.op.poisson_jacobi',	-- krylov is having errors.  TODO bug in its boundary code?
			}
			self.solver.ops:insert(NoDiv{
				solver = self.solver,
				potentialField = 'mPot',
				
				--[=[ using div (m/rho) = 0, solve for div m:
				vectorField = 'm',
			
				-- div v = 0
				-- div (m/ρ) = 0
				-- 1/ρ div m - 1/ρ^2 m dot grad ρ = 0
				-- div m = (m dot grad ρ)/ρ 
				chargeCode = self:template[[
	<? for j=0,solver.dim-1 do ?>{
		global <?=cons_t?> const * const Ujm = U - solver->stepsize.s<?=j?>;
		global <?=cons_t?> const * const Ujp = U + solver->stepsize.s<?=j?>;
		real const drho_dx = (Ujp->rho - Ujm->rho) * (.5 / solver->grid_dx.s<?=j?>);
		source -= drho_dx * U->m.s<?=j?> / U->rho;
	}<? end ?>
]],
				--]=]
				-- [=[ reading via div(v), writing via div(m)
				readVectorField = function(op,offset,j)
					local function U(field) return 'U['..offset..'].'..field end
					return U('m.s'..j)..' / '..U'rho'
				end,
				writeVectorField = function(op,dv)
					return self:template([[
#if 0	// just adjust velocity
	U->m = real3_sub(U->m, real3_real_mul(<?=dv?>, U->rho));
#endif

#if 0	// adjust ETotal as well
	U->ETotal -= .5 * U->rho * coordLenSq(U->m, pt);
	U->m = real3_sub(U->m, real3_real_mul(<?=dv?>, U->rho));
	U->ETotal += .5 * U->rho * coordLenSq(U->m, pt);
#endif

#if 1 	// recalculate cons
//// MODULE_DEPENDS: <?=primFromCons?> <?=consFromPrim?>
	<?=prim_t?> W;
	<?=primFromCons?>(&W, solver, U, pt);
	W.v = real3_sub(W.v, <?=dv?>);
	<?=consFromPrim?>(U, solver, &W, pt);
#endif
]], {dv=dv})
				end,
				codeDepends = {
					assert(self.symbols.eqn_common),
				},
				--]=]
			})
		end
	end
end

function MHD:getSymbolFields()
	return MHD.super.getSymbolFields(self):append{
		'roe_t',
		'calcRoeValues',
		'eigen_forRoeAvgs',
	}
end

--[[
B^2 = kg^2/(C^2 s^2)
B^2 / mu = kg^2/(C^2 s^2) * C^2/(kg*m) = kg/(m s^2)
--]]
MHD.guiVars = table{
	{name='heatCapacityRatio', value=2},	-- 5/3 for most problems, but 2 for Brio-Wu, so I will just set it here for now (in case something else is reading it before it is set there)
	
	-- works good as mu0 = 1
	-- solver->mu0 / unit_kg_m_per_C2 == 1
	-- solver->mu0 / (unit_kg * unit_m / (unit_C * unit_C)) == 1
	-- unit_C = sqrt((unit_kg * unit_m) / solver->mu0);
	{name='mu0', value=constants.vacuumPermeability_in_kg_m_per_C2, units='(kg*m)/C^2'},
	
	{name='Cp', value=1, compileTime=true},
}

if MHD.useFixedCh then
	MHD.guiVars:insert{name='Ch', value=.1}
end

function MHD:initCodeModules()
	MHD.super.initCodeModules(self)

	self.solver.modules:add{
		name = self.symbols.roe_t,
		structs = {self.roeStruct},
		-- only generated for cl, not for ffi cdef
		headercode = 'typedef '..self.symbols.roe_t..' roe_t;',
	}
	
	-- where do I put this to make it the default value for MHD solvers,
	-- but not override a value set by the init state?
	self.guiVars.coulomb.value = math.sqrt(self.guiVars.kilogram.value * self.guiVars.meter.value / self.guiVars.mu0.value)
end

-- don't use default
function MHD:initCodeModule_fluxFromCons() end
function MHD:initCodeModule_consFromPrim_primFromCons() end

MHD.solverCodeFile = 'hydro/eqn/glm-mhd.cl'

MHD.displayVarCodeUsesPrims = true

MHD.predefinedDisplayVars = {
	'U rho',
	'U m',
	'U ETotal',
	'U B',
	'U div B',
	'U psi',
}

function MHD:getDisplayVars()
	local vars = MHD.super.getDisplayVars(self)
	vars:append{
		{name='v', code='value.vreal3 = W.v;', type='real3', units='m/s'},
		{name='P', code='value.vreal = W.P;', units='kg/(m*s^2)'},
		{name='PMag', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_PMag(solver, &W, x);
]], units='kg/(m*s^2)'},
		{name='PTotal', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = W.P + calc_PMag(solver, &W, x);
]], units='kg/(m*s^2)'},
		{name='eInt', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_eInt(solver, &W);
]], units='m^2/s^2'},
		{name='EInt', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_EInt(solver, &W);
]], units='kg/(m*s^2)'},
		{name='eKin', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_eKin(&W, x);
]], units='m^2/s^2'},
		{name='EKin', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_EKin(&W, x);
]], units='kg/(m*s^2)'},
		{name='eHydro', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_eHydro(solver, &W, x);
]], units='m^2/s^2'},
		{name='EHydro', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_EHydro(solver, &W, x);
]], units='kg/(m*s^2)'},
		{name='EM energy', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_EM_energy(solver, &W, x);
]], units='kg/(m*s^2)'},
		{name='eTotal', code='value.vreal = U->ETotal / W.rho;', units='m^2/s^2'},
		{name='S', code='value.vreal = W.P / pow(W.rho, (real)solver->heatCapacityRatio);'},
		{name='H', code='value.vreal = calc_H(solver, W.P);', units='kg/(m*s^2)'},
		{name='h', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_H(solver, W.P) / W.rho;
]], units='m^2/s^2'},
		{name='HTotal', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_HTotal(solver, &W, U->ETotal, x);
]], units='kg/(m*s^2)'},
		{name='hTotal', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_hTotal(solver, &W, U->ETotal, x);
]], units='m^2/s^2'},
		{name='speed of sound', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = calc_Cs(solver, &W);
]], units='m/s'},
		{name='alfven velocity', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal3 = calc_CA(solver, U);
]], type='real3', units='m/s'},
		{name='Mach number', code = self:template[[
//// MODULE_DEPENDS: <?=eqn_common?>
value.vreal = coordLen(W.v, x) / calc_Cs(solver, &W);
]]},
		{name='temperature', code=self:template[[
<? local materials = require 'hydro.materials' ?>
#define C_v				<?=('%.50f'):format(materials.Air.C_v)?>
//// MODULE_DEPENDS: <?=eqn_common?>
	value.vreal = calc_eInt(solver, &W) / C_v;
]], units='K'},
		{name='primitive reconstruction error', code=self:template[[
//prim have just been reconstructed from cons
//so reconstruct cons from prims again and calculate the difference
<?=cons_t?> U2;
//// MODULE_DEPENDS: <?=consFromPrim?>
<?=consFromPrim?>(&U2, solver, &W, x);
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
if (!<?=OOB?>(1,1)) {
//// MODULE_DEPENDS: <?=eqn.gravOp.symbols.calcGravityAccel?>
	<?=eqn.gravOp.symbols.calcGravityAccel?>(&value.vreal3, solver, U, x);
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

function MHD:eigenWaveCode(args)
	local eig = '('..args.eig..')'
	return ({
		eig..'->v.x - '..eig..'->Cf',
		eig..'->v.x - '..eig..'->CAx',
		eig..'->v.x - '..eig..'->Cs',
		eig..'->v.x',
		eig..'->v.x + '..eig..'->Cs',
		eig..'->v.x + '..eig..'->CAx',
		eig..'->v.x + '..eig..'->Cf',
		
		--#warning there's a few PLM routines that expect eigenvalues to be ordered ... so replace them with a eigen_calcMinMaxWaves
		'-'..eig..'->Ch',
		eig..'->Ch',
	})[args.waveIndex+1] or error("got a bad waveIndex")
end

-- TODO consWaveCode[Prefix] for plm

function MHD:consWaveCodeMinMax(args)
	return self:template([[
range_t lambda;
<?=calcCellMinMaxEigenvalues?>(&lambda, solver, <?=U?>, <?=pt?>, <?=n?>); 
<?=eqn:waveCodeAssignMinMax(declare, resultMin, resultMax, 'lambda.min', 'lambda.max')?>
]], args)
end

-- since there is no prefix code
MHD.consWaveCodeMinMaxAllSides = MHD.consWaveCodeMinMax

return MHD

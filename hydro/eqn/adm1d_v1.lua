--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" 2008 chapter on Toy 1+1 spacetimes.

See comments in my gravitation-waves project adm1d_v1.lua file for the math.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local symmath = require 'symmath'
local EinsteinEqn = require 'hydro.eqn.einstein'

local ADM_BonaMasso_1D_2008Alcubierre = class(EinsteinEqn)

ADM_BonaMasso_1D_2008Alcubierre.name = 'ADM_BonaMasso_1D_2008Alcubierre' 

ADM_BonaMasso_1D_2008Alcubierre.consVars = {
	{name='alpha', type='real'}, 
	{name='gamma_xx', type='real'}, 
	{name='a_x', type='real'}, 
	{name='D_g', type='real'}, 
	{name='KTilde', type='real'},
}
ADM_BonaMasso_1D_2008Alcubierre.numWaves = 3	-- alpha and gamma_xx are source-term only

--ADM_BonaMasso_1D_2008Alcubierre.hasFluxFromConsCode = true
ADM_BonaMasso_1D_2008Alcubierre.useSourceTerm = true
ADM_BonaMasso_1D_2008Alcubierre.roeUseFluxFromCons = true


ADM_BonaMasso_1D_2008Alcubierre.guiVars = {
	{name='a_x_convCoeff', value=10},
	{name='D_g_convCoeff', value=10},
}

-- code that goes in initCond and in the solver
function ADM_BonaMasso_1D_2008Alcubierre:getCommonFuncCode()
	return template([[
void setFlatSpace(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* U,
	real3 x
) {
	*U = (<?=eqn.cons_t?>){
		.alpha = 1, 
		.gamma_xx = 1,
		.a_x = 0,
		.D_g = 0,
		.KTilde = 0,
	};
}

]], {
		eqn = self,
		solver = self.solver,
	})
end

ADM_BonaMasso_1D_2008Alcubierre.initCondCode = [[
kernel void applyInitCond(
	constant <?=solver.solver_t?>* solver,
	constant <?=solver.initCond_t?>* initCond,
	global <?=eqn.cons_t?>* UBuf,
	const global <?=coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	real3 x = cellBuf[index].pos;
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = sym3_ident;
	sym3 K_ll = sym3_zero;

	<?=code?>

	U->alpha = alpha;
	U->gamma_xx = gamma_ll.xx;
	U->KTilde = K_ll.xx / sqrt(gamma_ll.xx);
}

kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	real dx_alpha = (U[1].alpha - U[-1].alpha) / solver->grid_dx.x;
	real dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / solver->grid_dx.x;

	U->a_x = dx_alpha / U->alpha;
	U->D_g = dx_gamma_xx / U->gamma_xx;
}
]]

ADM_BonaMasso_1D_2008Alcubierre.solverCodeFile = 'hydro/eqn/adm1d_v1.cl'

function ADM_BonaMasso_1D_2008Alcubierre:getDisplayVars()
	return ADM_BonaMasso_1D_2008Alcubierre.super.getDisplayVars(self):append{
		-- adm1d_v2 cons vars:
		{name='d_xxx', code='value.vreal = .5 * U->D_g * U->gamma_xx;'},
		{name='K_xx', code='value.vreal = U->KTilde * sqrt(U->gamma_xx);'},
		-- aux:
		{name='dx_alpha', code='value.vreal = U->alpha * U->a_x;'},
		{name='dx_gamma_xx', code='value.vreal = U->gamma_xx * U->D_g;'},
		{name='volume', code='value.vreal = U->alpha * sqrt(U->gamma_xx);'},
		{name='f', code='value.vreal = calc_f(U->alpha);'},
		{name='df/dalpha', code='value.vreal = calc_dalpha_f(U->alpha);'},
		{name='K', code='value.vreal = U->KTilde / sqrt(U->gamma_xx);'},
		{name='expansion', code='value.vreal = -U->KTilde / sqrt(U->gamma_xx);'},
		{name='gravity mag', code='value.vreal = -U->alpha * U->alpha * U->a_x / U->gamma_xx;'},
	
		{name='alpha vs a_x', code=[[
	if (OOB(1,1)) {
		value.vreal = 0.;
	} else {
		real dx_alpha = (U[1].alpha - U[-1].alpha) / (2. * solver->grid_dx.x);
		value.vreal = fabs(dx_alpha - U->alpha * U->a_x);
	}
]]},

		{name='gamma_xx vs D_g', code=[[
	if (OOB(1,1)) {
		value.vreal = 0.;
	} else {
		real dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / (2. * solver->grid_dx.x);
		value.vreal = fabs(dx_gamma_xx - U->gamma_xx * U->D_g);
	}
]]},
	}
end

ADM_BonaMasso_1D_2008Alcubierre.eigenVars = table{
	{name='f', type='real'},
	{name='alpha', type='real'},
	{name='gamma_xx', type='real'},
}

function ADM_BonaMasso_1D_2008Alcubierre:eigenWaveCodePrefix(n, eig, x, waveIndex)
	return template([[
	real eig_lambda = <?=eig?>.alpha * sqrt(<?=eig?>.f / <?=eig?>.gamma_xx);
]], {
		eig = '('..eig..')',
	})
end

function ADM_BonaMasso_1D_2008Alcubierre:eigenWaveCode(n, eig, x, waveIndex)
	if waveIndex == 0 then
		return '-eig_lambda'
	elseif waveIndex == 1 then
		return '0'
	elseif waveIndex == 2 then
		return 'eig_lambda'
	else
		error'got a bad waveIndex'
	end
end

function ADM_BonaMasso_1D_2008Alcubierre:consWaveCodePrefix(n, U, x, waveIndex)
	return template([[
	real f = calc_f(<?=U?>.alpha);
	real eig_lambda = <?=U?>.alpha * sqrt(f / <?=U?>.gamma_xx);
]], {
		U = '('..U..')',
	})
end

ADM_BonaMasso_1D_2008Alcubierre.consWaveCode = ADM_BonaMasso_1D_2008Alcubierre.eigenWaveCode
	

-- TODO store flat values somewhere, then perturb all real values here
--  then you can move this into the parent class
local function crand() return 2 * math.random() - 1 end
function ADM_BonaMasso_1D_2008Alcubierre:fillRandom(epsilon)
	local ptr = ADM_BonaMasso_1D_2008Alcubierre.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.numCells-1 do
		ptr[i].alpha = ptr[i].alpha + 1
		ptr[i].gamma_xx = ptr[i].gamma_xx + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return ADM_BonaMasso_1D_2008Alcubierre

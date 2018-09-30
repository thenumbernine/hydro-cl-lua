--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" 2008 chapter on Toy 1+1 spacetimes.

See comments in my gravitation-waves project adm1d_v1.lua file for the math.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local symmath = require 'symmath'
local EinsteinEqn = require 'eqn.einstein'

local ADM_BonaMasso_1D_2008Alcubierre = class(EinsteinEqn)

ADM_BonaMasso_1D_2008Alcubierre.name = 'ADM_BonaMasso_1D_2008Alcubierre' 

ADM_BonaMasso_1D_2008Alcubierre.consVars = {
	{alpha = 'real'}, 
	{gamma_xx = 'real'}, 
	{a_x = 'real'}, 
	{D_g = 'real'}, 
	{KTilde = 'real'},
}
ADM_BonaMasso_1D_2008Alcubierre.numWaves = 3	-- alpha and gamma_xx are source-term only

ADM_BonaMasso_1D_2008Alcubierre.mirrorVars = {{'gamma_xx', 'a_x', 'D_g', 'KTilde'}}

ADM_BonaMasso_1D_2008Alcubierre.hasEigenCode = true
--ADM_BonaMasso_1D_2008Alcubierre.hasFluxFromConsCode = true
ADM_BonaMasso_1D_2008Alcubierre.useSourceTerm = true
ADM_BonaMasso_1D_2008Alcubierre.roeUseFluxFromCons = true


ADM_BonaMasso_1D_2008Alcubierre.guiVars = {
	{name='a_x_convCoeff', value=10},
	{name='D_g_convCoeff', value=10},
}

-- code that goes in initState and in the solver
function ADM_BonaMasso_1D_2008Alcubierre:getCommonFuncCode()
	return template([[
void setFlatSpace(global <?=eqn.cons_t?>* U, real3 x) {
	*U = (<?=eqn.cons_t?>){
		.alpha = 1, 
		.gamma_xx = 1,
		.a_x = 0,
		.D_g = 0,
		.KTilde = 0,
	};
}
]], {eqn=self})
end

ADM_BonaMasso_1D_2008Alcubierre.initStateCode = [[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(mins, maxs), .5);
	
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

ADM_BonaMasso_1D_2008Alcubierre.solverCodeFile = 'eqn/adm1d_v1.cl'

function ADM_BonaMasso_1D_2008Alcubierre:getDisplayVars()
	return ADM_BonaMasso_1D_2008Alcubierre.super.getDisplayVars(self):append{
		-- adm1d_v2 cons vars:
		{d_xxx = '*value = .5 * U->D_g * U->gamma_xx;'},
		{K_xx = '*value = U->KTilde * sqrt(U->gamma_xx);'},
		-- aux:
		{dx_alpha = '*value = U->alpha * U->a_x;'},
		{dx_gamma_xx = '*value = U->gamma_xx * U->D_g;'},
		{volume = '*value = U->alpha * sqrt(U->gamma_xx);'},
		{f = '*value = calc_f(U->alpha);'},
		{['df/dalpha'] = '*value = calc_dalpha_f(U->alpha);'},
		{K = '*value = U->KTilde / sqrt(U->gamma_xx);'},
		{expansion = '*value = -U->KTilde / sqrt(U->gamma_xx);'},
		{['gravity mag'] = '*value = -U->alpha * U->alpha * U->a_x / U->gamma_xx;'},
	
		{['alpha vs a_x'] = [[
	if (OOB(1,1)) {
		*value = 0.;
	} else {
		real dx_alpha = (U[1].alpha - U[-1].alpha) / (2. * solver->grid_dx.x);
		*value = fabs(dx_alpha - U->alpha * U->a_x);
	}
]]},

		{['gamma_xx vs D_g'] = [[
	if (OOB(1,1)) {
		*value = 0.;
	} else {
		real dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / (2. * solver->grid_dx.x);
		*value = fabs(dx_gamma_xx - U->gamma_xx * U->D_g);
	}
]]},
	}
end

ADM_BonaMasso_1D_2008Alcubierre.eigenVars = table{
	{f = 'real'},
	{alpha = 'real'},
	{gamma_xx = 'real'},
}

function ADM_BonaMasso_1D_2008Alcubierre:eigenWaveCodePrefix(side, eig, x, waveIndex)
	return template([[
	real eig_lambda = <?=eig?>.alpha * sqrt(<?=eig?>.f / <?=eig?>.gamma_xx);
]], {
		eig = '('..eig..')',
	})
end

function ADM_BonaMasso_1D_2008Alcubierre:eigenWaveCode(side, eig, x, waveIndex)
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

function ADM_BonaMasso_1D_2008Alcubierre:consWaveCodePrefix(side, U, x, waveIndex)
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

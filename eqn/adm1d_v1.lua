--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" 2008 chapter on Toy 1+1 spacetimes.

See comments in my gravitation-waves project adm1d_v1.lua file for the math.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local symmath = require 'symmath'
local NumRelEqn = require 'eqn.numrel'

local ADM_BonaMasso_1D_Alcubierre2008 = class(NumRelEqn)

ADM_BonaMasso_1D_Alcubierre2008.name = 'ADM_BonaMasso_1D_Alcubierre2008' 

ADM_BonaMasso_1D_Alcubierre2008.consVars = {'alpha', 'gamma_xx', 'a_x', 'D_g', 'KTilde'}
ADM_BonaMasso_1D_Alcubierre2008.numWaves = 3	-- alpha and gamma_xx are source-term only

ADM_BonaMasso_1D_Alcubierre2008.mirrorVars = {{'gamma_xx', 'a_x', 'D_g', 'KTilde'}}

ADM_BonaMasso_1D_Alcubierre2008.hasEigenCode = true
ADM_BonaMasso_1D_Alcubierre2008.useSourceTerm = true
ADM_BonaMasso_1D_Alcubierre2008.hasFluxFromCons = true

function ADM_BonaMasso_1D_Alcubierre2008:getCodePrefix()
	-- pick out whatever variables that 'codes' needs to convert
	local lines = table()

		-- don't call super because it generates the guivar code
		-- which is already being generated in initState
		--ADM_BonaMasso_1D_Alcubierre2008.super.getCodePrefix(self),
		
	lines:insert(template([[
void setFlatSpace(global <?=eqn.cons_t?>* U) {
	*U = (<?=eqn.cons_t?>){
		.alpha = 1, 
		.gamma_xx = 1,
		.a_x = 0,
		.D_g = 0,
		.KTilde = 0,
	};
}
]], {eqn=self}))
		
	if self.initState.getCodePrefix then
		lines:insert(self.initState:getCodePrefix(self.solver, function(exprs, vars)
			return {
				alpha  = exprs.alpha,
				gamma_xx = exprs.gamma[1],	-- only need g_xx
				a_x = (exprs.alpha:diff(vars[1]) / exprs.alpha)(),	-- only need a_x
				D_g = (exprs.gamma[1]:diff(vars[1]) / exprs.gamma[1])(),	-- only need D_xxx
				KTilde = exprs.K[1] / symmath.sqrt(exprs.gamma[1]),	-- only need K_xx
			}
		end))
	end

	return lines:concat'\n'
end

ADM_BonaMasso_1D_Alcubierre2008.initStateCode = [[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	U->alpha = calc_alpha(x.x, x.y, x.z);
	U->gamma_xx = calc_gamma_xx(x.x, x.y, x.z);
	U->a_x = calc_a_x(x.x, x.y, x.z);
	U->D_g = calc_D_g(x.x, x.y, x.z);
	U->KTilde = calc_KTilde(x.x, x.y, x.z);
}
]]

function ADM_BonaMasso_1D_Alcubierre2008:getSolverCode()
	return template(file['eqn/adm1d_v1.cl'], {eqn=self, solver=self.solver})
end

function ADM_BonaMasso_1D_Alcubierre2008:getDisplayVars()
	return {
		-- source-only:
		{alpha = '*value = U->alpha;'},
		{gamma_xx = '*value = U->gamma_xx;'},
		-- both 1998 and 2008 cons vars:
		{a_x = '*value = U->a_x;'},
		-- 1998-only cons vars:
		{d_xxx = '*value = .5 * U->D_g * U->gamma_xx;'},
		{K_xx = '*value = U->KTilde * sqrt(U->gamma_xx);'},
		-- 2008-only cons vars:
		{D_g = '*value = U->D_g;'},
		{KTilde = '*value = U->KTilde;'},
		-- aux:
		{dx_alpha = '*value = U->alpha * U->a_x;'},
		{dx_gamma_xx = '*value = U->gamma_xx * U->D_g;'},
		{volume = '*value = U->alpha * sqrt(U->gamma_xx);'},
		{f = '*value = calc_f(U->alpha);'},
		{['df/dalpha'] = '*value = calc_dalpha_f(U->alpha);'},
	}
end

function ADM_BonaMasso_1D_Alcubierre2008:getEigenTypeCode()
	return template([[
typedef struct {
	real f, alpha, gamma_xx;
} <?=eqn.eigen_t?>;
]], {eqn=self, solver=self.solver})
end
		
function ADM_BonaMasso_1D_Alcubierre2008:getEigenDisplayVars()
	return {
		{f = '*value = eigen->f;'},
		{alpha = '*value = eigen->alpha;'},
		{gamma_xx = '*value = eigen->gamma_xx;'},
	}
end

local ffi = require 'ffi'
local function crand() return 2 * math.random() - 1 end
function ADM_BonaMasso_1D_Alcubierre2008:fillRandom(epsilon)
	local solver = self.solver
	local ptr = ffi.new(self.cons_t..'[?]', solver.volume)
	for i=0,solver.volume-1 do
		ptr[i].alpha = epsilon * crand()
		ptr[i].gamma_xx = 1 + epsilon * crand()
		ptr[i].a_x = epsilon * crand()
		ptr[i].D_g = epsilon * crand()
		ptr[i].KTilde = epsilon * crand()
	end
	solver.UBufObj:fromCPU(ptr)
end

return ADM_BonaMasso_1D_Alcubierre2008

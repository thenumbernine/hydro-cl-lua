--[[
Based on the Alcubierre 1997 "The appearance of coorindate shocks in hyperbolic formalisms of General Relativity".
This is also a 1D version of the 3D formalism in the Alcubierre 2008 book "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 

a_x,t + (α f K_xx / γ_xx),x = 0
d_xxx,t + (α K_xx),x = 0
K_xx,t + (α a_x),x = (α / γ_xx) (a_x d_xxx - K_xx^2)

for d_xxx = 1/2 γ_xx,x

expanded:

a_x,t + α,x f K_xx / γ_xx + α f,x K_xx / γ_xx + α f K_xx,x / γ_xx - α f K_xx / γ_xx^2 γ_xx,x = 0
a_x,t + α f / γ_xx K_xx,x = α K_xx / γ_xx (f (2 d_xxx / γ_xx - a_x) - α a_x f')

d_xxx,t + α,x K_xx + α K_xx,x = 0
d_xxx,t + α K_xx,x = -α a_x K_xx

K_xx,t + α,x a_x + α a_x,x = α / γ_xx (a_x d_xxx - K_xx^2)
K_xx,t + α a_x,x = α ((a_x d_xxx - K_xx^2) / γ_xx - a_x^2)

[ a_x ]     [ 0, 0, α f / γ_xx ] [ a_x ]     [ α K_xx / γ_xx (f (2 d_xxx / γ_xx - a_x) - α a_x f') ]
[d_xxx]   + [ 0, 0,      α     ] [d_xxx]   = [ -α a_x K_xx                                         ]
[ K_xx],t   [ α, 0,      0     ] [ K_xx],x   [ α ((a_x d_xxx - K_xx^2) / γ_xx - a_x^2)             ]

... has eigenvalues ...

Lambda = {-α √(f/γ_xx), 0, α √(f/γ_xx)}
   
... and eigenvectors ...

	[ √(f/γ_xx), 0, √(f/γ_xx) ]
Q = [ √(γ_xx/f), 1, √(γ_xx/f) ]
    [     -1,    0,     1     ]

       [ √(γ_xx/f)/2, 0, -1/2 ]
Q^-1 = [     -γ_xx/f, 1,   0  ]
       [ √(γ_xx/f)/2, 0,  1/2 ]

corresponding eigenfields:

( should be a_x - f d_xxx / γ_xx, √(f) K_xx / γ_xx \pm a_x / √(γ_xx) )

eigenfields:

[ √(γ_xx/f)/2, 0, -1/2 ] [ a_x,a ]
[     -γ_xx/f, 1,   0  ] [d_xxx,a]
[ √(γ_xx/f)/2, 0,  1/2 ] [ K_xx,a]

1/2 (√(γ_xx / f) a_x,a - K_xx,a) ... * -2 * √(f) / γ_xx
-γ_xx / f a_x,a + d_xxx,a ... * -f / γ_xx
1/2 (√(γ_xx / f) a_x,a + K_xx,a) ... * 2 * √(f) / γ_xx

√(f) / γ_xx K_xx,a - 1 / √(γ_xx) a_x,a	<- check
a_x,a - f / γ_xx * d_xxx,a 				<- check
√(f) / γ_xx K_xx,a + 1 / √(γ_xx) a_x,a	<- check




same system, favoring flux terms, incorporating α and γ to do just that ...
(I'm trying to figure out why adding in the extra source terms that come from linearizing wrt the primitive variables messes the equation up, but removing them works fine)

[  α  ]     [              0,                       0,          0,   0,      0     ] [  α  ]     [           -α^2 f K            ]
[γ_xx ]     [              0,                       0,          0,   0,      0     ] [γ_xx ]     [          -2 α K_xx            ]
[ a_x ]  +  [ α^2 K_xx / γ_xx (f + α f'), -α K_xx f / γ_xx^2,   0,   0, α f / γ_xx ] [ a_x ]   = [               0               ]
[d_xxx]     [          α^2 K_xx,                    0,          0,   0,      α     ] [d_xxx]     [               0               ]
[K_xx ],t   [           α^2 a_x,                    0,          α,   0,      0     ] [K_xx ],x   [ α / γ_xx (a_x d_xxx - K_xx^2) ]

...and voila, our source term now matches up with what the paper says.
the catch?  finding the eigenvectors is difficult, since the eigenvalues are ±α √(f/γ_xx) and 0 x3
for this reason I use the eigenfields, and it reconstructs the above matrix 
 all except the α and γ columns ...

so why can't those terms go into the source?
why do I have to use a matrix that reconstructs without them?

--]]

local class = require 'ext.class'
local table = require 'ext.table'
local EinsteinEqn = require 'hydro.eqn.einstein'

local ADM_BonaMasso_1D_1997Alcubierre = class(EinsteinEqn)
ADM_BonaMasso_1D_1997Alcubierre.name = 'ADM_BonaMasso_1D_1997_Alcubierre'

ADM_BonaMasso_1D_1997Alcubierre.numStates = 5
ADM_BonaMasso_1D_1997Alcubierre.numWaves = 3

ADM_BonaMasso_1D_1997Alcubierre.consVars = {
	{name='alpha', type='real'},
	{name='gamma_xx', type='real', variance=''},
	{name='a_x', type='real', variance=''}, 
	{name='d_xxx', type='real', variance=''}, 
	{name='K_xx', type='real', variance=''},
}

function ADM_BonaMasso_1D_1997Alcubierre:createInitState()
	ADM_BonaMasso_1D_1997Alcubierre.super.createInitState(self)
	self:addGuiVars{
		-- hmm, it is useful to make this proportional to dt ... since it's used for a constraint ...
		{name='a_x_convCoeff', value=10},
		{name='d_xxx_convCoeff', value=10},
	}
end

-- don't use default
function ADM_BonaMasso_1D_1997Alcubierre:initCodeModule_fluxFromCons() end

-- don't use eqn.einstein:
function ADM_BonaMasso_1D_1997Alcubierre:createDisplayComponents() end

ADM_BonaMasso_1D_1997Alcubierre.solverCodeFile = 'hydro/eqn/adm1d_v2.cl'

function ADM_BonaMasso_1D_1997Alcubierre:getDisplayVars()
	return ADM_BonaMasso_1D_1997Alcubierre.super.getDisplayVars(self):append{
		-- adm1d_v1 cons vars:
		{name='D_g', code='value.vreal = 2. * U->d_xxx / U->gamma_xx;'},
		{name='KTilde', code='value.vreal = U->K_xx / sqrt(U->gamma_xx);'},
		-- aux:
		{name='dx_alpha', code='value.vreal = U->alpha * U->a_x;'},
		{name='dx_gamma_xx', code='value.vreal = 2. * U->d_xxx;'},
		{name='volume', code='value.vreal = U->alpha * sqrt(U->gamma_xx);'},
		{name='f', code='value.vreal = calc_f(U->alpha);'},
		{name='df/dalpha', code='value.vreal = calc_dalpha_f(U->alpha);'},
		{name='K', code='value.vreal = U->K_xx / U->gamma_xx;'},
		{name='expansion', code='value.vreal = -U->K_xx / U->gamma_xx;'},
		{name='gravity mag', code='value.vreal = -U->alpha * U->alpha * U->a_x / U->gamma_xx;'},
	
		{name='alpha vs a_x', code=self:template[[
	if (<?=OOB?>(1,1)) {
		value.vreal = 0.;
	} else {
		real dx_alpha = (U[1].alpha - U[-1].alpha) / (2. * solver->grid_dx.x);
		value.vreal = fabs(dx_alpha - U->alpha * U->a_x);
	}
]]},

		{name='gamma_xx vs d_xxx', code=self:template[[
	if (<?=OOB?>(1,1)) {
		value.vreal = 0.;
	} else {
		real dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / (2. * solver->grid_dx.x);
		value.vreal = fabs(dx_gamma_xx - 2. * U->d_xxx);
	}
]]},
	}
end

ADM_BonaMasso_1D_1997Alcubierre.eigenVars = {
	{name='alpha', type='real'},
	{name='sqrt_f_over_gamma_xx', type='real'},
}

function ADM_BonaMasso_1D_1997Alcubierre:eigenWaveCodePrefix(args)
	return self:template([[
real const eig_lambda = (<?=eig?>)->alpha * (<?=eig?>)->sqrt_f_over_gamma_xx;
]], args)
end

function ADM_BonaMasso_1D_1997Alcubierre:eigenWaveCode(args)
	if args.waveIndex == 0 then
		return '-eig_lambda'
	elseif args.waveIndex == 1 then
		return '0'
	elseif args.waveIndex == 2 then
		return 'eig_lambda'
	else
		error'got a bad waveIndex'
	end
end

function ADM_BonaMasso_1D_1997Alcubierre:consWaveCodePrefix(args)
	return self:template([[
real const f = calc_f((<?=U?>)->alpha);
real const eig_lambda = (<?=U?>)->alpha * sqrt(f / (<?=U?>)->gamma_xx);
]], args)
end
ADM_BonaMasso_1D_1997Alcubierre.consWaveCode = ADM_BonaMasso_1D_1997Alcubierre.eigenWaveCode

return ADM_BonaMasso_1D_1997Alcubierre

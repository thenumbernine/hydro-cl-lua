--[[
Based on the Alcubierre 1997 "The appearance of coorindate shocks in hyperbolic formalisms of General Relativity".
This is also a 1D version of the 3D formalism in the Alcubierre 2008 book "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 

a_x,t + (alpha f K_xx / gamma_xx),x = 0
d_xxx,t + (alpha K_xx),x = 0
K_xx,t + (alpha a_x),x = (alpha / gamma_xx) (a_x d_xxx - K_xx^2)

for d_xxx = 1/2 gamma_xx,x

expanded:

a_x,t + alpha,x f K_xx / gamma_xx + alpha f,x K_xx / gamma_xx + alpha f K_xx,x / gamma_xx - alpha f K_xx / gamma_xx^2 gamma_xx,x = 0
a_x,t + alpha f / gamma_xx K_xx,x = alpha K_xx / gamma_xx (f (2 d_xxx / gamma_xx - a_x) - alpha a_x f')

d_xxx,t + alpha,x K_xx + alpha K_xx,x = 0
d_xxx,t + alpha K_xx,x = -alpha a_x K_xx

K_xx,t + alpha,x a_x + alpha a_x,x = alpha / gamma_xx (a_x d_xxx - K_xx^2)
K_xx,t + alpha a_x,x = alpha ((a_x d_xxx - K_xx^2) / gamma_xx - a_x^2)

[ a_x ]     [   0,   0, alpha f / gamma_xx ] [ a_x ]     [ alpha K_xx / gamma_xx (f (2 d_xxx / gamma_xx - a_x) - alpha a_x f') ]
[d_xxx]   + [   0,   0,        alpha       ] [d_xxx]   = [ -alpha a_x K_xx                                                     ]
[ K_xx],t   [ alpha, 0,          0         ] [ K_xx],x   [ alpha ((a_x d_xxx - K_xx^2) / gamma_xx - a_x^2)                     ]

... has eigenvalues ...

Lambda = {-alpha sqrt(f/gamma_xx), 0, alpha sqrt(f/gamma_xx)}
   
... and eigenvectors ...

	[ sqrt(f/gamma_xx), 0, sqrt(f/gamma_xx) ]
Q = [ sqrt(gamma_xx/f), 1, sqrt(gamma_xx/f) ]
    [       -1,         0,         1        ]

       [ sqrt(gamma_xx/f)/2, 0, -1/2 ]
Q^-1 = [     -gamma_xx/f,    1,   0  ]
       [ sqrt(gamma_xx/f)/2, 0,  1/2 ]

corresponding eigenfields:

( should be a_x - f d_xxx / gamma_xx, sqrt(f) K_xx / gamma_xx \pm a_x / sqrt(gamma_xx) )

eigenfields:

[ sqrt(gamma_xx/f)/2, 0, -1/2 ] [ a_x,a ]
[     -gamma_xx/f,    1,   0  ] [d_xxx,a]
[ sqrt(gamma_xx/f)/2, 0,  1/2 ] [ K_xx,a]

1/2 (sqrt(gamma_xx / f) a_x,a - K_xx,a) ... * -2 * sqrt(f) / gamma_xx
-gamma_xx / f a_x,a + d_xxx,a ... * -f / gamma_xx
1/2 (sqrt(gamma_xx / f) a_x,a + K_xx,a) ... * 2 * sqrt(f) / gamma_xx

sqrt(f) / gamma_xx K_xx,a - 1 / sqrt(gamma_xx) a_x,a	<- check
a_x,a - f / gamma_xx * d_xxx,a 							<- check
sqrt(f) / gamma_xx K_xx,a + 1 / sqrt(gamma_xx) a_x,a	<- check




same system, favoring flux terms, incorporating alpha and gamma to do just that ...
(I'm trying to figure out why adding in the extra source terms that come from linearizing wrt the primitive variables messes the equation up, but removing them works fine)

[  alpha ]     [                      0,                               0,              0,   0,          0         ] [  alpha ]     [              -alpha^2 f K             ]
[gamma_xx]     [                      0,                               0,              0,   0,          0         ] [gamma_xx]     [             -2 alpha K_xx             ]
[   a_x  ]  +  [ alpha^2 K_xx / gamma_xx (f + alpha f'), -alpha K_xx f / gamma_xx^2,   0,   0, alpha f / gamma_xx ] [   a_x  ]   = [                    0                  ]
[  d_xxx ]     [                 alpha^2 K_xx,                         0,              0,   0,        alpha       ] [  d_xxx ]     [                    0                  ]
[  K_xx  ],t   [                  alpha^2 a_x,                         0,            alpha, 0,          0         ] [  K_xx  ],x   [ alpha / gamma_xx (a_x d_xxx - K_xx^2) ]

...and voila, our source term now matches up with what the paper says.
the catch?  finding the eigenvectors is difficult, since the eigenvalues are +-alpha sqrt(f/gamma_xx) and 0 x3
for this reason I use the eigenfields, and it reconstructs the above matrix 
 all except the alpha and gamma columns ...

so why can't those terms go into the source?
why do I have to use a matrix that reconstructs without them?

--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local template = require 'template'
local EinsteinEqn = require 'eqn.einstein'

local ADM_BonaMasso_1D_Alcubierre1997 = class(EinsteinEqn)
ADM_BonaMasso_1D_Alcubierre1997.name = 'ADM_BonaMasso_1D_Alcubierre1997'

ADM_BonaMasso_1D_Alcubierre1997.numStates = 5
ADM_BonaMasso_1D_Alcubierre1997.numWaves = 3

ADM_BonaMasso_1D_Alcubierre1997.consVars = {
	{alpha = 'real'},
	{gamma_xx = 'real'},
	{a_x = 'real'}, 
	{d_xxx = 'real'}, 
	{K_xx = 'real'},
}
ADM_BonaMasso_1D_Alcubierre1997.mirrorVars = {{'gamma_xx', 'a_x', 'd_xxx', 'K_xx'}}

ADM_BonaMasso_1D_Alcubierre1997.hasEigenCode = true
ADM_BonaMasso_1D_Alcubierre1997.useSourceTerm = true

function ADM_BonaMasso_1D_Alcubierre1997:createInitState()
	ADM_BonaMasso_1D_Alcubierre1997.super.createInitState(self)
	self:addGuiVars{
		-- hmm, it is useful to make this proportional to dt ... since it's used for a constraint ...
		{name='a_x_convCoeff', value=10},
		{name='d_xxx_convCoeff', value=10},
	}
end

function ADM_BonaMasso_1D_Alcubierre1997:getCommonFuncCode()
	return template([[
void setFlatSpace(global <?=eqn.cons_t?>* U) {
	*U = (<?=eqn.cons_t?>){
		.alpha = 1, 
		.gamma_xx = 1,
		.a_x = 0,
		.d_xxx = 0,
		.K_xx = 0,
	};
}
]], {eqn=self})
end

ADM_BonaMasso_1D_Alcubierre1997.initStateCode = [[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real alpha = 1.;
	real3 beta_u = _real3(0,0,0);
	sym3 gamma_ll = _sym3(1,0,0,1,0,1);
	sym3 K_ll = _sym3(0,0,0,0,0,0);

	<?=code?>

	U->alpha = alpha;
	U->gamma_xx = gamma_ll.xx;
	U->K_xx = K_ll.xx;
}

kernel void initDerivs(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	real dx_alpha = (U[1].alpha - U[-1].alpha) / grid_dx0;
	real dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / grid_dx0;
	
	U->a_x = dx_alpha / U->alpha;
	U->d_xxx = .5 * dx_gamma_xx;
}
]]

function ADM_BonaMasso_1D_Alcubierre1997:getSolverCode()
	return template(file['eqn/adm1d_v2.cl'], {eqn=self, solver=self.solver})
end

function ADM_BonaMasso_1D_Alcubierre1997:getDisplayVars()
	return ADM_BonaMasso_1D_Alcubierre1997.super.getDisplayVars(self):append{
		-- adm1d_v1 cons vars:
		{D_g = '*value = 2. * U->d_xxx / U->gamma_xx;'},
		{KTilde = '*value = U->K_xx / sqrt(U->gamma_xx);'},
		-- aux:
		{dx_alpha = '*value = U->alpha * U->a_x;'},
		{dx_gamma_xx = '*value = 2. * U->d_xxx;'},
		{volume = '*value = U->alpha * sqrt(U->gamma_xx);'},
		{f = '*value = calc_f(U->alpha);'},
		{['df/dalpha'] = '*value = calc_dalpha_f(U->alpha);'},
		{K = '*value = U->K_xx / U->gamma_xx;'},
		{expansion = '*value = -U->K_xx / U->gamma_xx;'},
		{['gravity mag'] = '*value = -U->alpha * U->alpha * U->a_x / U->gamma_xx;'},
	
		{['alpha vs a_x'] = [[
	if (OOB(1,1)) {
		*value = 0.;
	} else {
		real dx_alpha = (U[1].alpha - U[-1].alpha) / (2. * grid_dx0);
		*value = fabs(dx_alpha - U->alpha * U->a_x);
	}
]]},

		{['gamma_xx vs d_xxx'] = [[
	if (OOB(1,1)) {
		*value = 0.;
	} else {
		real dx_gamma_xx = (U[1].gamma_xx - U[-1].gamma_xx) / (2. * grid_dx0);
		*value = fabs(dx_gamma_xx - 2. * U->d_xxx);
	}
]]},
	}
end

ADM_BonaMasso_1D_Alcubierre1997.eigenVars = {
	{alpha = 'real'},
	{sqrt_f_over_gamma_xx = 'real'},
}

function ADM_BonaMasso_1D_Alcubierre1997:eigenWaveCodePrefix(side, eig, x, waveIndex)
	return template([[
	real eig_lambda = <?=eig?>->alpha * <?=eig?>->sqrt_f_over_gamma_xx;
]], {
		eig = '('..eig..')',
	})
end

function ADM_BonaMasso_1D_Alcubierre1997:eigenWaveCode(side, eig, x, waveIndex)
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

function ADM_BonaMasso_1D_Alcubierre1997:fillRandom(epsilon)
	local ptr = ADM_BonaMasso_1D_Alcubierre1997.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.volume-1 do
		ptr[i].alpha = ptr[i].alpha + 1
		ptr[i].gamma_xx = ptr[i].gamma_xx + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return ADM_BonaMasso_1D_Alcubierre1997

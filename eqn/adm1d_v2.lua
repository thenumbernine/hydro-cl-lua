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

--]]

local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'eqn.eqn'

local ADM_BonaMasso_1D_Alcubierre1997 = class(Equation)
ADM_BonaMasso_1D_Alcubierre1997.name = 'ADM_BonaMasso_1D_Alcubierre1997'

ADM_BonaMasso_1D_Alcubierre1997.numStates = 5
ADM_BonaMasso_1D_Alcubierre1997.numWaves = 3

ADM_BonaMasso_1D_Alcubierre1997.consVars = {'alpha', 'gamma_xx', 'a_x', 'd_xxx', 'K_xx'}
ADM_BonaMasso_1D_Alcubierre1997.mirrorVars = {{'gamma_xx', 'a_x', 'd_xxx', 'K_xx'}}

ADM_BonaMasso_1D_Alcubierre1997.useSourceTerm = true

ADM_BonaMasso_1D_Alcubierre1997.displayVars = table()
	:append(ADM_BonaMasso_1D_Alcubierre1997.consVars)
	:append{'dx_alpha', 'dx_gamma_xx', 'D_g', 'KTilde_xx', 'volume'}

ADM_BonaMasso_1D_Alcubierre1997.initStates = require 'eqn.init_adm'
ADM_BonaMasso_1D_Alcubierre1997.initStateNames = table.map(ADM_BonaMasso_1D_Alcubierre1997.initStates, function(state) return state.name end)

function ADM_BonaMasso_1D_Alcubierre1997:codePrefix()
	return table.map(self.codes, function(code,name,t)
		return 'real calc_'..name..code, #t+1
	end):concat'\n'
end

ADM_BonaMasso_1D_Alcubierre1997.guiVars = {'f'}
ADM_BonaMasso_1D_Alcubierre1997.f = {
	value = 0,	-- 0-based index into options
	name = 'f',
	options = {'1', '1.69', '.49', '1 + 1/alpha^2'},
}

function ADM_BonaMasso_1D_Alcubierre1997:getInitStateCode(solver)
	local initState = self.initStates[solver.initStatePtr[0]+1]
	
	local alphaVar = require 'symmath'.var'alpha'
	self.codes = initState.init(solver, ({
		{f = 1},
		{f = 1.69},
		{f = 1.49},
		{f = 1 + 1/alphaVar^2, alphaVar=alphaVar},
	})[self.f.value+1])

	return table{
		self:codePrefix(),
		[[
__kernel void initState(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real4 x = CELL_X(i);
	__global cons_t* U = UBuf + index;
	
	U->alpha = calc_alpha(x.x, x.y, x.z);
	U->gamma_xx = calc_gamma_xx(x.x, x.y, x.z);
	U->a_x = calc_a_x(x.x, x.y, x.z);
	U->d_xxx = calc_d_xxx(x.x, x.y, x.z);
	U->K_xx = calc_K_xx(x.x, x.y, x.z);
}
]],
	}:concat'\n'
end

function ADM_BonaMasso_1D_Alcubierre1997:solverCode()
	return table{
		self:codePrefix(),
		'#include "eqn/adm1d_v2.cl"',
	}:concat'\n'
end

ADM_BonaMasso_1D_Alcubierre1997.eigenVars = {'sqrt_g_xx_over_f'}
function ADM_BonaMasso_1D_Alcubierre1997:getEigenInfo()
	local makeStruct = require 'eqn.makestruct'
	return {
		typeCode =
			makeStruct('eigen_t', self.eigenVars) .. '\n' ..
			makeStruct('fluxXform_t', {'alpha', 'gamma_xx', 'f'}),
		code = nil,
		displayVars = self.eigenVars,
	}
end

return ADM_BonaMasso_1D_Alcubierre1997 

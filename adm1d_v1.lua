--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" 2008 chapter on Toy 1+1 spacetimes.

conservative variables:
a_x = (ln alpha),x = alpha,x / alpha
D_g = (ln gamma_xx),x = gamma_xx,x / gamma_xx
KTilde_xx = sqrt(gamma_xx) K_xx

d_xxx = 1/2 gamma_xx,x
D_g = 2 d_xxx / gamma_xx
d_xxx = 1/2 D_g gamma_xx

K_xx = KTilde_xx / sqrt(gamma_xx)
KTilde_xx,x = K_xx,x sqrt(gamma_xx) + 1/2 sqrt(gamma_xx) K_xx gamma_xx,x / gamma_xx
	= K_xx,x sqrt(gamma_xx) + 1/2 KTilde_xx D_g
K_xx,x = KTilde_xx,x / sqrt(gamma_xx) - 1/2 KTilde_xx D_g / sqrt(gamma_xx)

flux form:
a_x,t + (alpha f K_xx),x = 0
D_g,t + (2 alpha K_xx),x = 0
KTilde_xx,t + (alpha a_x / sqrt(gamma_xx)),x = 0

expanded derivatives:
a_x,t + alpha,x f K_xx + alpha f,x K_xx + alpha f K_xx,x = 0
a_x,t + alpha f / sqrt(gamma_xx) KTilde_xx,x = alpha KTilde_xx / sqrt(gamma_xx) (f (1/2 D_g - a_x) - a_x alpha f'))

D_g,t + 2 alpha,x K_xx + 2 alpha K_xx,x = 0
D_g,t + 2 alpha / sqrt(gamma_xx) KTilde_xx,x = 2 alpha KTilde_xx / sqrt(gamma_xx) (1/2 D_g - a_x)

KTilde_xx,t + alpha,x a_x / sqrt(gamma_xx) + alpha a_x,x / sqrt(gamma_xx) + alpha a_x * -1/2 1/sqrt(gamma_xx)^3 gamma_xx,x = 0
KTilde_xx,t + alpha / sqrt(gamma_xx) a_x,x = alpha a_x / sqrt(gamma_xx) (1/2 D_g - a_x)

[   a_x   ]     [         0,              0, alpha f / sqrt(gamma_xx) ] [a_x ]     [ alpha KTilde_xx / sqrt(gamma_xx) (f (1/2 D_g - a_x) - a_x alpha f') ]
[   D_g   ]   + [         0,              0, 2 alpha / sqrt(gamma_xx) ] [D_g ]   = [        2 alpha KTilde_xx / sqrt(gamma_xx) (1/2 D_g - a_x)           ]
[KTilde_xx],t   [ alpha / sqrt(gamma_xx), 0,            0             ] [K_xx],x   [           alpha a_x / sqrt(gamma_xx) (1/2 D_g - a_x)                ]

... has eigenvalues ...

Lambda = {-alpha sqrt(f/gamma_xx), 0, alpha sqrt(f/gamma_xx)}
   
... and eigenvectors ...

	[     f,    0,    f    ]
Q = [     2,    1,    2    ]
    [ -sqrt(f), 0, sqrt(f) ]

       [ 1/(2f), 0, -1/(2 sqrt(f)) ]
Q^-1 = [ -2/f,   1,        0       ]
       [ 1/(2f), 0,  1/(2 sqrt(f)) ]

--]]

local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local ADM_BonaMasso_1D_Alcubierre2008 = class(Equation)
ADM_BonaMasso_1D_Alcubierre2008.name = 'ADM_BonaMasso_1D_Alcubierre2008' 

ADM_BonaMasso_1D_Alcubierre2008.numStates = 5
ADM_BonaMasso_1D_Alcubierre2008.numWaves = 3	-- alpha and gamma_xx are source-term only

ADM_BonaMasso_1D_Alcubierre2008.consVars = {'alpha', 'gamma_xx', 'a_x', 'D_g', 'KTilde_xx'}
ADM_BonaMasso_1D_Alcubierre2008.mirrorVars = {{'gamma_xx', 'a_x', 'D_g', 'KTilde_xx'}}

ADM_BonaMasso_1D_Alcubierre2008.useSourceTerm = true

ADM_BonaMasso_1D_Alcubierre2008.displayVars = table()
	:append(ADM_BonaMasso_1D_Alcubierre2008.consVars)
	:append{'dx_alpha', 'dx_gamma_xx', 'd_xxx', 'K_xx', 'volume'}

ADM_BonaMasso_1D_Alcubierre2008.initStates = require 'init_adm'
ADM_BonaMasso_1D_Alcubierre2008.initStateNames = table.map(ADM_BonaMasso_1D_Alcubierre2008.initStates, function(state) return state.name end)

function ADM_BonaMasso_1D_Alcubierre2008:codePrefix()
	return table.map(self.codes, function(code,name,t)
		return 'real calc_'..name..code, #t+1
	end):concat'\n'
end

ADM_BonaMasso_1D_Alcubierre2008.guiVars = {'f'}
ADM_BonaMasso_1D_Alcubierre2008.f = {
	value = 0,	-- 0-based index into options
	name = 'f',
	options = {'1', '1.69', '.49', '1 + 1/alpha^2'},
}

function ADM_BonaMasso_1D_Alcubierre2008:getInitStateCode(solver)
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
	U->D_g = (2. * calc_d_xxx(x.x, x.y, x.z) / U->gamma_xx);
	U->KTilde_xx = calc_K_xx(x.x, x.y, x.z) * sqrt(U->gamma_xx);
}
]],
	}:concat'\n'
end

function ADM_BonaMasso_1D_Alcubierre2008:solverCode()
	return table{
		self:codePrefix(),
		'#include "adm1d_v1.cl"',
	}:concat'\n'
end

ADM_BonaMasso_1D_Alcubierre2008.eigenVars = {'f'}
function ADM_BonaMasso_1D_Alcubierre2008:getEigenInfo()
	local makeStruct = require 'makestruct'
	return {
		typeCode =
			makeStruct('eigen_t', self.eigenVars) .. '\n' ..
			makeStruct('fluxXform_t', {'alpha', 'gamma_xx', 'f'}),
		code = nil,
		displayVars = self.eigenVars,
	}
end

return ADM_BonaMasso_1D_Alcubierre2008 

local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local ADM1D3to5Var = class(Equation)
ADM1D3to5Var.name = 'ADM1D3to5Var' 

ADM1D3to5Var.numStates = 5

ADM1D3to5Var.consVars = {'alpha', 'gamma_xx', 'a_x', 'd_xxx', 'KTilde_xx'}
ADM1D3to5Var.initStates = {'gaussianWave'}

ADM1D3to5Var.displayVars = table()
	:append(ADM1D3to5Var.consVars)
	:append{'K_xx', 'volume'}

function ADM1D3to5Var:solverCode(clnumber)
	local symmath = require 'symmath'
symmath.tostring = require 'symmath.tostring.SingleLine'		
	local x = symmath.var'x'
	local alphaVar = symmath.var'alpha'
	local xc = 150
	local H = 5
	local sigma = 10
	local h = H * symmath.exp(-(x - xc)^2 / sigma^2)
	
	local alpha = 1
	local gamma_xx = 1 - h:diff(x)^2
	local K_xx = -h:diff(x,x) / gamma_xx^.5

	local kappa = 1
	local f = 1 + kappa / alphaVar^2
	
	local function makesym(expr, k)
		return symmath.clone(expr):simplify(), k
	end
	
	local exprs = table{
		alpha = alpha,
		gamma_xx = gamma_xx,
		K_xx = K_xx,
	}:map(makesym)
	exprs.dx_alpha = exprs.alpha:diff(x):simplify()
	exprs.dx_gamma_xx = exprs.gamma_xx:diff(x):simplify()

	local function compileC(expr, args)
print('compile',expr,table.map(args,tostring):concat', ')
		local code = require 'symmath.tostring.Lua'(expr, args)
print('...to',code)
		return code
	end

	local calc = exprs:map(function(expr, name)
		return compileC(expr, {x}), name
	end)
	
	f = makesym(f)
	calc.f = compileC(f, {alphaVar})
	
	local dalpha_f = f:diff(alphaVar):simplify()
	calc.dalpha_f = compileC(dalpha_f, {alphaVar})

	return '#include "adm1d3to5var.cl"'
end

function ADM1D3to5Var:getEigenTypeCode()
	return require 'makestruct'(self.eigenType, {'f'})
end
-- override / don't use the default matrix stuff
function ADM1D3to5Var:getEigenCode() end

return ADM1D3to5Var 

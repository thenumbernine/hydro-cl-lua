--[[
From 2010 Colliander, Simpso, Sulem - Numerical Simulations of the Energy-Supercritical nonlinear Schrodinger equation
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'eqn.eqn'

local NLSEqn = class(Equation)
NLSEqn.name = 'NLSEqn'

NLSEqn.hasCalcDTCode = true

-- TODO just use cplx_t?
NLSEqn.consVars = {{re='real'}, {im='real'}}

NLSEqn.initStates = require 'init.nls'

NLSEqn.initStateCode = [[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);

	real r = x.x;
	real re = 0;
	real im = 0;
	
	<?=code?>

	UBuf[index] = (<?=eqn.cons_t?>){.re=re, .im=im};
}
]]

NLSEqn.solverCodeFile = 'eqn/nls.cl'

function NLSEqn:getDisplayVars()
	return NLSEqn.super.getDisplayVars(self):append{
		{norm = '*value = sqrt(U->re*U->re + U->im*U->im);'},
	}
end

return NLSEqn 

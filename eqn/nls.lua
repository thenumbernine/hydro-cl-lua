local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'eqn.eqn'

local NLSEqn = class(Equation)
NLSEqn.name = 'NLSEqn'
NLSEqn.consVars = {'re', 'im'}

NLSEqn.initStates = require 'init.nls'

NLSEqn.initStateCode = [[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);

	real r = x.x;
	real re = 0;
	real im = 0;
	
	<?=code?>

	UBuf[index] = (<?=eqn.cons_t?>){.re=re, .im=im};
}
]]

NLSEqn.solverCodeFile = 'eqn/nls.cl'

function NLSEqn:getDisplayVars()
	return {
		{re = '*value = U->re;'},
		{im = '*value = U->im;'},
		{norm = '*value = sqrt(U->re*U->re + U->im*U->im);'},
	}
end

return NLSEqn 

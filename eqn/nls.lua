--[[
From 2010 Colliander, Simpso, Sulem - Numerical Simulations of the Energy-Supercritical nonlinear Schrodinger equation
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'eqn.eqn'

local NLSEqn = class(Equation)
NLSEqn.name = 'NLSEqn'

NLSEqn.hasCalcDTCode = true

-- this tells eqn to not provide its own eigen code
NLSEqn.hasEigenCode = true
NLSEqn.hasFluxFromConsCode = true

-- TODO just use cplx_t?
NLSEqn.consVars = {{re='real'}, {im='real'}}

function NLSEqn:createBoundaryOptions()
	self.solver.boundaryOptions:insert{
		oscillating = function(args)
			local index = args.index
			local assign = args.assign
			local lines = table()
			local gridSizeSide = 'gridSize_'..xNames[args.side]
			if args.minmax == 'min' then
				-- buf[0] = buf[4]
				-- buf[1] = buf[3]
				-- buf[2] is ... ?
				lines:insert(
					assign('buf['..index'j'..']',
						'buf['..index('2*numGhost-j')..']'))
			elseif args.minmax == 'max' then
				lines:insert(
					assign('buf['..index(gridSizeSide..'-1-j')..']',
						'buf['..index(gridSizeSide..'-2*numGhost-2+j')..']'))
			end
			return lines:concat'\n'
		end,
	}
	self.solver.boundaryOptions:insert{
		fixed = function(args)
			if args.minmax == 'min' then
				assign('buf['..index'j'..']', '0')
			elseif args.minmax == 'max' then
				assign('buf['..index(gridSizeSide'-1-j')..']', '0')
			end
		end,
	}
end

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

-- SolverBase adds eqn:getEigenTypeCode(), even though it's predominantly a finite volume / Godunov (Roe) solver property
function NLSEqn:getEigenTypeCode() end

return NLSEqn 

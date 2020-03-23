--[[
From 2010 Colliander, Simpso, Sulem - Numerical Simulations of the Energy-Supercritical nonlinear Schrodinger equation
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'eqn.eqn'

local common = require 'common'
local xNames = common.xNames

local NLSEqn = class(Equation)
NLSEqn.name = 'nls'

NLSEqn.hasCalcDTCode = true

-- this tells eqn to not provide its own eigen code
NLSEqn.hasEigenCode = true
NLSEqn.hasFluxFromConsCode = true
NLSEqn.useSourceTerm = true

-- TODO just use cplx_t?
NLSEqn.consVars = {
	{name='q', type='cplx'},
}

function NLSEqn:createBoundaryOptions()
	local Boundary = self.solver.Boudary
	local BoundaryOscillating = class(Boundary)
	BoundaryOscillating.name = 'oscillating'
	function BoundaryOscillating:getCode(args)
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
	end
	self.solver:addBoundaryOption(BoundaryOscillating)
	
	local BoundaryFixed = class(Boundary)	
	BoundaryFixed.name = 'fixed'
	function BoundaryFixed:getCode(args)
		if args.minmax == 'min' then
			assign('buf['..index'j'..']', '0')
		elseif args.minmax == 'max' then
			assign('buf['..index(gridSizeSide'-1-j')..']', '0')
		end
	end
	self.solver:addBoundaryOption(BoundaryFixed)
end

NLSEqn.predefinedDisplayVars = {
	'U q re',
	'U q im',
	'U q abs',
	'U q arg',
}

NLSEqn.initStates = require 'init.nls'

NLSEqn.initStateCode = [[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);

	real r = fabs(x.x);
	real re = 0;
	real im = 0;
	
	<?=code?>

	UBuf[index].q = _cplx(re, im);
}
]]

NLSEqn.solverCodeFile = 'eqn/nls.cl'

-- SolverBase adds eqn:getEigenTypeCode(), and eqn/eqn.lua provides a default, so this overrides it 
function NLSEqn:getEigenTypeCode() end

return NLSEqn 

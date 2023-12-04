--[[
From 2010 Colliander, Simpso, Sulem - Numerical Simulations of the Energy-Supercritical nonlinear Schrodinger equation
--]]
local table = require 'ext.table'
local Equation = require 'hydro.eqn.eqn'

local xNames = require 'hydro.common'.xNames

local NLSEqn = Equation:subclass()
NLSEqn.name = 'nls'

-- TODO just use cplx_t?
NLSEqn.consVars = {
	{name='q', type='cplx'},
}

NLSEqn.solverCodeFile = 'hydro/eqn/nls.cl'
NLSEqn.initConds = require 'hydro.init.nls':getList()

function NLSEqn:createBoundaryOptions()
	local Boundary = self.solver.Boudary
	local BoundaryOscillating = Boundary:subclass()
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
					'buf['..index('2 * solver->numGhost - j')..']'))
		elseif args.minmax == 'max' then
			lines:insert(
				assign('buf['..index(gridSizeSide..'-1-j')..']',
					'buf['..index(gridSizeSide..' - 2 * solver->numGhost - 2 + j')..']'))
		end
		return lines:concat'\n'
	end
	self.solver:addBoundaryOption(BoundaryOscillating)
	
	local BoundaryFixed = Boundary:subclass()	
	BoundaryFixed.name = 'fixed'
	function BoundaryFixed:getCode(args)
		local index = args.index
		local assign = args.assign
		if args.minmax == 'min' then
			assign('buf['..index'j'..']', '0')
		elseif args.minmax == 'max' then
			assign('buf['..index(gridSizeSide'-1-j')..']', '0')
		end
	end
	self.solver:addBoundaryOption(BoundaryFixed)
end

-- the default display-all is broken since i switched to the pick-component option
NLSEqn.predefinedDisplayVars = {
	'U q re',
	'U q im',
	'U q abs',
	'U q arg',
}

-- don't use default
function NLSEqn:initCodeModule_calcDT() end

return NLSEqn 

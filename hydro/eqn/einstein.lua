--[[
common functions for all Einstein field equation solvers
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local Equation = require 'hydro.eqn.eqn'

local common = require 'hydro.common'
local xNames = common.xNames

local EinsteinEquation = class(Equation)

EinsteinEquation.initConds = table():append(
	require 'hydro.init.senr':getList()
):append(
	require 'hydro.init.einstein':getList()
)

function EinsteinEquation:createInitState()
	EinsteinEquation.super.createInitState(self)
	self:addGuiVars{
		{
			type = 'combo',
			name = 'f_eqn',
			options = {
				'2/alpha',							-- 1+log slicing
				'1 + 1/alpha^2', 					-- Alcubierre 10.2.24: "shock avoiding condition" for Toy 1+1 spacetimes 
				'1', 								-- Alcubierre 4.2.50 - harmonic slicing
				'0', '.49', '.5', '1.5', '1.69',
			}
		},
	}
end

function EinsteinEquation:initCodeModules()
	EinsteinEquation.super.initCodeModules(self)

	self.solver.modules:addFromMarkup(self:template(file['hydro/eqn/einstein.cl']))
end

-- add an option for fixed Minkowsky boundary spacetime
-- TODO now there is already a BoundaryFixed in hydro/solver/gridsolver, but no easy way to parameterize how to set what fixed values it is
function EinsteinEquation:createBoundaryOptions()
	local Boundary = self.solver.Boudary
	local BoundaryFixed = class(Boundary)
	BoundaryFixed.name = 'fixed'
	function BoundaryFixed:getCode(args)
		local lines = table()
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		for _,j in ipairs{'j', gridSizeSide..'-numGhost+j'} do
			local index = args.indexv(j)
			local U = 'buf[INDEX('..index..')]'
			lines:insert(self:template([[
	setFlatSpace(solver, &<?=U?>, cell_x((int4)(<?=index?>, 0)));
]], 		{
				U = U,
				index = index,
			}))
		end
		return lines:concat'\n'
	end
	
	self.solver:addBoundaryOption(BoundaryFixed)
end

function EinsteinEquation:getModuleDepends_displayCode() 
	return {
		-- for the addDisplayComponents 
		'calc_gamma_ll',
		'calc_gamma_uu',
	}
end

-- any modules this code needs, add to solver module dependencies
function EinsteinEquation:createDisplayComponents()
	local solver = self.solver
	solver:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gamma_ij',
		code = [[
const global <?=eqn.cons_t?>* U = buf + index;
sym3 gamma_ll = calc_gamma_ll(U, x);
value->vreal = real3_weightedLen(value->vreal3, gamma_ll);
]],
	})
	solver:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gamma^ij',
		code = [[
const global <?=eqn.cons_t?>* U = buf + index;
sym3 gamma_uu = calc_gamma_uu(U, x);
value->vreal = real3_weightedLen(value->vreal3, gamma_uu);
]],
	})
	solver:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted gamma^ij',
		code = [[
const global <?=eqn.cons_t?>* U = buf + index;
sym3 gamma_uu = calc_gamma_uu(U, x);
value->vreal = sym3_dot(value->vsym3, gamma_uu);]],
	})
end

return EinsteinEquation

--[[
common functions for all Einstein field equation solvers
--]]

local table = require 'ext.table'
local path = require 'ext.path'
local Equation = require 'hydro.eqn.eqn'

local common = require 'hydro.common'
local xNames = common.xNames

local EinsteinEquation = Equation:subclass()

EinsteinEquation.initConds = table(require 'hydro.init.einstein':getList())
do
	local success, initConds = pcall(require, 'hydro.init.senr')
	if success then
		EinsteinEquation.initConds:append(initConds:getList())
	end
end

function EinsteinEquation:getSymbolFields() 
	return EinsteinEquation.super.getSymbolFields(self):append{
		-- defined in subclasses:
		'setFlatSpace',
		'calc_gamma_ll',
		'calc_gamma_uu',

		-- placeholder for rescale[From/To]Coord_*
		'rescaleFromCoord_rescaleToCoord',
		'cplx3_rescaleFromCoord_cplx3_rescaleToCoord',
		'_3sym3_rescaleFromCoord__3sym3_rescaleToCoord',
		'sym3sym3_rescaleFromCoord_sym3sym3_rescaleToCoord',
	}
end

function EinsteinEquation:createInitState()
	EinsteinEquation.super.createInitState(self)
	self:addGuiVars{
		{
			type = 'combo',
			name = 'f_eqn',
			options = {
				'2/alpha',							-- 1+log slicing
				'1 + 1/alpha^2', 					-- 2008 Alcubierre 10.2.24: "shock avoiding condition" for Toy 1+1 spacetimes 
				'1', 								-- 2008 Alcubierre 4.2.50 - harmonic slicing
				'0', '.49', '.5', '1.5', '1.69',
			}
		},
	}
end

function EinsteinEquation:initCodeModules()
	EinsteinEquation.super.initCodeModules(self)

	self.solver.modules:addFromMarkup(self:template(path'hydro/eqn/einstein.cl':read()))
end

-- add an option for fixed Minkowsky boundary spacetime
-- TODO now there is already a BoundaryFixed in hydro/solver/gridsolver, but no easy way to parameterize how to set what fixed values it is
function EinsteinEquation:createBoundaryOptions()
	local eqn = self
	local Boundary = self.solver.Boudary
	local BoundaryFixed = Boundary:subclass()
	BoundaryFixed.name = 'fixed'
	function BoundaryFixed:getCode(args)
		local lines = table()
		local gridSizeSide = 'solver->gridSize.'..xNames[args.side]
		for _,j in ipairs{'j', gridSizeSide..' - solver->numGhost + j'} do
			local index = args.indexv(j)
			lines:insert(eqn:template([[
<?=setFlatSpace?>(solver, &buf[INDEX(<?=index?>)], cellBuf[<?=index?>].pos);
]], 		{
				index = index,
			}))
		end
		return lines:concat'\n', {eqn.symbols.setFlatSpace}
	end
	
	self.solver:addBoundaryOption(BoundaryFixed)
end

-- convenient trackvars:
-- TODO generate these by expression in CL automatically
function EinsteinEquation:getDisplayVars()
	local vars = EinsteinEquation.super.getDisplayVars(self)
	
	vars:append{
		{name='alpha-1', code=[[ value.vreal = U->alpha - 1.;]], type='real'},
		{
			name = 'gamma_ll',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_gamma_ll?>
value.vsym3 = <?=calc_gamma_ll?>(U, x);
]],
			type = 'sym3',
		},
		{
			name = 'gamma_uu',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_gamma_uu?>
value.vsym3 = <?=calc_gamma_uu?>(U, x);
]],
			type = 'sym3',
		},
	}

	return vars
end

-- any modules this code needs, add to solver module dependencies
function EinsteinEquation:createDisplayComponents()
	local solver = self.solver
	solver:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gamma_ij',
		code = self:template[[
global <?=cons_t?> const * const U = buf + index;
sym3 gamma_ll = <?=calc_gamma_ll?>(U, x);
value->vreal = real3_weightedLen(value->vreal3, gamma_ll);
]],
	})
	solver:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gamma^ij',
		code = self:template[[
global <?=cons_t?> const * const U = buf + index;
sym3 gamma_uu = <?=calc_gamma_uu?>(U, x);
value->vreal = real3_weightedLen(value->vreal3, gamma_uu);
]],
	})
	solver:addDisplayComponent('sym3', {
		onlyFor = 'U',
		name = 'tr weighted gamma^ij',
		code = self:template[[
global <?=cons_t?> const * const U = buf + index;
sym3 gamma_uu = <?=calc_gamma_uu?>(U, x);
value->vreal = sym3_dot(value->vsym3, gamma_uu);]],
	})
end

return EinsteinEquation

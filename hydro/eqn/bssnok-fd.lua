local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local EinsteinEqn = require 'hydro.eqn.einstein'

local common = require 'hydro.common'
local xNames = common.xNames


local BSSNOKFiniteDifferenceEquationBase = class(EinsteinEqn)

-- seems all the hyperbolic formalisms listed in Alcubierre's book use alpha sqrt(gamma^ii) for the speed-of-light wavespeed
-- however the 2017 Ruchlin paper says to use gamma_ij
--BSSNOKFiniteDifferenceEquationBase.cflMethod = '2008 Alcubierre'
--BSSNOKFiniteDifferenceEquationBase.cflMethod = '2013 Baumgarte et al, eqn 32'
BSSNOKFiniteDifferenceEquationBase.cflMethod = '2017 Ruchlin et al, eqn 53'

function BSSNOKFiniteDifferenceEquationBase:init(args)
	self.cflMethod = args.cflMethod
	BSSNOKFiniteDifferenceEquationBase.super.init(self, args)
end

function BSSNOKFiniteDifferenceEquationBase:getSymbolFields()
	return BSSNOKFiniteDifferenceEquationBase.super.getSymbolFields(self):append{
		'calc_det_gammaHat',
		'calc_det_gammaBar',
		'calc_det_gammaBarLL',
		'calc_gammaHat_ll',
		'calc_gammaHat_uu',
		'calc_gammaHat_LL',
		'calc_gammaHat_UU',
		'calc_gammaBar_ll',
		'calc_gammaBar_LL',
		'calc_gammaBar_uu',
		'calc_gammaBar_UU',
		'calc_exp_neg4phi',
		'mystery_C_U',

		-- placeholder for solver/bssnok-fd-pirk
		'BSSNOK_PIRK',
	}
end

-- source: https://en.wikipedia.org/wiki/Finite_difference_coefficient 
-- derivCoeffs[derivative][order] = {coeffs...}
local derivCoeffs = {
	-- upwind 1st deriv coefficients
	{
		[2] = {[0] = -1, 1},
		[4] = {[0] = -3, 4, -1, denom=2},
		[6] = {[0] = -11, 18, -9, 2, denom=6},
--[[		
		[8] = {[0] = -25, 48, -36, 16, -3, denom=12},
		senr_8 = { 	-- SENR's coeffs
			[-1] = -3,
			[0] = -10,
			[1] = 18,
			[2] = -6,
			[3] = 1,
			denom = 12,
		},
--]]	
-- [[ using senr's 8 by default
		[8] = { 	-- SENR's coeffs
			[-1] = -3,
			[0] = -10,
			[1] = 18,
			[2] = -6,
			[3] = 1,
			denom = 12,
		},
		senr_8 = { 	-- SENR's coeffs
			[-1] = -3,
			[0] = -10,
			[1] = 18,
			[2] = -6,
			[3] = 1,
			denom = 12,
		},
--]]
	}
}

-- used by the BSSNOKFiniteDifference solvers
-- maybe a parent class for just them?
-- TODO higher orders plz
-- computes the partial between U[stepSize] and U[0]
function BSSNOKFiniteDifferenceEquationBase:makePartialUpwind(field, fieldType, nameOverride)
	local clnumber = require 'cl.obj.number'
	fieldType = fieldType or self:fieldTypeForVar(field)

	local suffix = 'l'
	if not field:find'_' then suffix = '_' .. suffix end

	local add, real_mul, zero
	if fieldType == 'real' then
		--[[ no codegen precedence
		add = function(x,y) return '(('..x..') + ('..y..'))' end
		real_mul = function(x,y) return '(('..x..') * ('..y..'))' end
		--]]
		-- [[ codegen with precedence
		add = function(x,y,xt,yt)
			return x..' + '..y, 'add'
		end
		real_mul = function(x,y,xt,yt)
			if xt == 'add' then x = '('..x..')' end
			if yt == 'add' then y = '('..y..')' end
			return x..' * '..y, 'mul'
		end
		--]]
		zero = '0.'
	else
		add = function(x,y) return fieldType..'_add('..x..', '..y..')' end
		real_mul = function(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
		zero = fieldType..'_zero'
	end

	local name = nameOverride or ('partial_'..field..suffix..'_upwind')
	local deriv = 1
	local order = 'senr_8'
	local d1coeffs = assert(derivCoeffs[deriv][order], "couldn't find d/dx^"..(deriv/2).." coefficients of order "..order)

	local lines = table()
	if fieldType == 'real' then
		lines:insert('\treal3 '..name..';')
	elseif fieldType == 'cplx' then
		lines:insert('\tcplx3 '..name..';')
	elseif fieldType == 'real3' then
		lines:insert('\treal3x3 '..name..';')
	elseif fieldType == 'cplx3' then
		lines:insert('\tcplx3x3 '..name..';')
	else
		lines:insert('\t'..fieldType..' '..name..'[3];')
	end
	for i,xi in ipairs(xNames) do
		local namei
		if fieldType == 'real' 
		or fieldType == 'cplx'
		or fieldType == 'real3'
		or fieldType == 'cplx3'
		then
			namei = name..'.'..xi
		else
			namei = name..'['..(i-1)..']'
		end
		local expr, exprtype = zero, 'value'
		if i <= self.solver.dim then
			local U = 'U->'..field
			expr, exprtype = real_mul(U, clnumber(d1coeffs[0]), 'value', 'value')
			for j,coeff in pairs(d1coeffs) do
			--for j=1,#d1coeffs do
				--local coeff = d1coeffs[j]
				if type(j) == 'number' and j ~= 0 then
					local U = 'U['..j..' * updir.'..xi..' * solver->stepsize.'..xi..'].'..field
				
					local subexpr, subexprtype = real_mul(U, clnumber(coeff))
					expr, exprtype = add(expr, subexpr, exprtype, subexprtype)
				end
			end
			local denom = 'solver->grid_dx.'..xi
			if d1coeffs.denom then
				denom = denom .. ' * ' .. clnumber(d1coeffs.denom)
			end
			expr, exprtype = real_mul(expr, '(1. / ('..denom..'))', exprtype, 'value')
		end
		expr, exprtype = real_mul(expr, '(real)updir.'..xi, exprtype, 'value')
		lines:insert('\t'..namei..' = '..expr..';')
	end
	return lines:concat'\n'
end

function BSSNOKFiniteDifferenceEquationBase:initCodeModules()
	BSSNOKFiniteDifferenceEquationBase.super.initCodeModules(self)

	self.solver.modules:addFromMarkup(self:template(file'hydro/eqn/bssnok-fd.cl':read()))
end

function BSSNOKFiniteDifferenceEquationBase:initCodeModule_calcDTCell()
	self.solver.modules:add{
		name = self.symbols.calcDTCell,
		depends = table{
			self.symbols.eqn_common,
		}:append(
			({
				['2008 Alcubierre'] = { 
					self.symbols.calc_gamma_uu,
				},
				['2017 Ruchlin et al, eqn 53'] = {
					self.solver.coord.symbols.coord_sqrt_gHol_ll_ij,
				},
			})[self.cflMethod]
		),
		code = self:template[[
void <?=calcDTCell?>(
	global real * const dt,
	constant <?=solver_t?> const * const solver,
	global <?=cons_t?> const * const U,
	global <?=cell_t?> const * const cell
) {
	real3 const x = cell->pos;

<? if eqn.cflMethod == '2008 Alcubierre' then
?>	real3s3 gamma_uu = <?=calc_gamma_uu?>(U, x);
<? end 
?>
	<? for side=0,solver.dim-1 do ?>{
<? 
if eqn.cflMethod == '2013 Baumgarte et al, eqn 32' then
	-- TODO if the grid is static then this only needs to be done once
	if side == 0 then 
?>		*(dt) = (real)min(*(dt), solver->grid_dx.x);
<?	elseif side == 1 then 
?>		*(dt) = (real)min(*(dt), .5 * solver->grid_dx.x * solver->grid_dx.y);
<? 	elseif side == 2 then 
?>		*(dt) = (real)min(*(dt), .5 * solver->grid_dx.x * sin(.5 * solver->grid_dx.y) * solver->grid_dx.z);
<? 	end 
else
	if eqn.cflMethod == '2008 Alcubierre' then 
?>		//this is asserting alpha and W >0, which they should be
		real absLambdaMax = U->alpha * sqrt(gamma_uu.<?=sym(side+1,side+1)?>);
		real dx = solver->grid_dx.s<?=side?>;
		*(dt) = (real)min(*(dt), dx / absLambdaMax);
<? 	elseif eqn.cflMethod == '2017 Ruchlin et al, eqn 53' then 
?>		// for wavespeeds alpha sqrt(gammaBar^ii)
		// and if we assume alpha > 1
		// and gamma_ii ~ 1/gamma^ii 
		// and gammaBar_ij ~ gammaHat_ij 
		// then the typical CFL equation: dt <= dx / lambda, lambda = alpha sqrt(gammaBar^ii)
		// turns into the SENR code: dt <= sqrt(gammaHat_ii) * dx
		real sqrt_gammaHat_ii = coord_sqrt_gHol_ll<?=side..side?>(x);
		real ds = sqrt_gammaHat_ii  * solver->grid_dx.s<?=side?>;
		*(dt) = (real)min(*(dt), ds);
<? 	else
		error("unknown cflMethod "..tostring(eqn.cflMethod))
	end
end
?>	}<? end ?>
}
]],
	}
end

function BSSNOKFiniteDifferenceEquationBase:getDisplayVars()
	local vars = BSSNOKFiniteDifferenceEquationBase.super.getDisplayVars(self)

	vars:append{
		{
			name = 'f',
			code = self:template[[
//// MODULE_DEPENDS: <?=initCond_codeprefix?>
value.vreal = calc_f(U->alpha);
]],
		},
		{
			name = 'df/dalpha',
			code = self:template[[
//// MODULE_DEPENDS: <?=initCond_codeprefix?>
value.vreal = calc_dalpha_f(U->alpha);
]],
		},
		
		{name='W-1', code=[[value.vreal = U->W - 1.;]], type='real'},
		{name='alpha-W', code=[[value.vreal = U->alpha - U->W;]], type='real'},		-- this is the "pre-collapse" initial condition for SENR UIUC
		
		{
			name = 'gammaHat_ll',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_gammaHat_ll?>
value.vreal3s3 = <?=calc_gammaHat_ll?>(x);
]],
			type = 'real3s3',
		},
		{
			name = 'gammaHat_uu',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_gammaHat_uu?>
value.vreal3s3 = <?=calc_gammaHat_uu?>(x);
]],
			type = 'real3s3',
		},
		{
			name = 'gammaBar_ll',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_gammaBar_ll?>
value.vreal3s3 = <?=calc_gammaBar_ll?>(U, x);
]],
			type = 'real3s3',
		},
		{
			name = 'gammaBar_uu',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_gammaBar_uu?>
value.vreal3s3 = <?=calc_gammaBar_uu?>(U, x);
]],
			type = 'real3s3',
		},
		{
			name = 'gammaBar_LL',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_gammaBar_LL?>
value.vreal3s3 = <?=calc_gammaBar_LL?>(U, x);
]],
			type = 'real3s3',
		},
		{
			name = 'gammaBar_UU',
			code = self:template[[
//// MODULE_DEPENDS: <?=calc_gammaBar_UU?>
value.vreal3s3 = <?=calc_gammaBar_UU?>(U, x);
]],
			type = 'real3s3',
		},
		{
			name = 'K_ll',
			code = self:template[[
real exp_4phi = 1. / <?=calc_exp_neg4phi?>(U);
real3s3 gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);
value.vreal3s3 = real3s3_real_mul(
	real3s3_add(
		real3s3_rescaleToCoord_LL(U->ABar_LL, x),
		real3s3_real_mul(gammaBar_ll, U->K / 3.)
	), exp_4phi);
]],
			type = 'real3s3',
		},

		{name='det gammaBar - det gammaHat', code = self:template[[
value.vreal = determinant(<?=calc_gammaBar_ll?>(U, x)) - <?=calc_det_gammaBar?>(x);
]]},
		{name='det gamma based on phi', code = self:template[[
real exp_neg4phi = <?=calc_exp_neg4phi?>(U);
real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
//// MODULE_DEPENDS: <?=calc_det_gammaHat?>	
real det_gamma = exp_12phi * <?=calc_det_gammaHat?>(x);
value.vreal = det_gamma;
]]},
	}

	return vars
end

function BSSNOKFiniteDifferenceEquationBase:createDisplayComponents()
	BSSNOKFiniteDifferenceEquationBase.super.createDisplayComponents(self)
	
	local solver = self.solver
	solver:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gammaBar_IJ',
		code = self:template[[
int index = INDEXV(solver, i);
real3 x = cellBuf[index].pos;
global <?=cons_t?> const * const U = buf + index;
real3s3 gammaBar_LL = <?=calc_gammaBar_LL?>(U, x);
value->vreal = real3_weightedLen(value->vreal3, gammaBar_LL);]],
	})
	solver:addDisplayComponent('real3', {
		onlyFor = 'U',
		name = 'norm weighted gammaBar_ij',
		code = self:template[[
int index = INDEXV(solver, i);
real3 x = cellBuf[index].pos;
global <?=cons_t?> const * const U = buf + index;
real3s3 gammaBar_ll = <?=calc_gammaBar_ll?>(U, x);
value->vreal = real3_weightedLen(value->vreal3, gammaBar_ll);]],
	})
	solver:addDisplayComponent('real3s3', {
		onlyFor = 'U',
		name = 'tr weighted gamma^IJ',
		code = self:template[[
int index = INDEXV(solver, i);
real3 x = cellBuf[index].pos;
global <?=cons_t?> const * const U = buf + index;
real3s3 gamma_UU = real3s3_rescaleFromCoord_uu(<?=calc_gamma_uu?>(U, x), x);
value->vreal = real3s3_dot(value->vreal3s3, gamma_UU);]],
	})
	solver:addDisplayComponent('real3s3', {
		onlyFor = 'U',
		name = 'tr weighted gammaBar^IJ',
		code = self:template[[
int index = INDEXV(solver, i);
real3 x = cellBuf[index].pos;
global <?=cons_t?> const * const U = buf + index;
real3s3 gammaBar_UU = <?=calc_gammaBar_UU?>(U, x);
value->vreal = real3s3_dot(value->vreal3s3, gammaBar_UU);]],
	})
end

return BSSNOKFiniteDifferenceEquationBase 

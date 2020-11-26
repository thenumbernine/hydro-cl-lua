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

	self.solver.modules:addFromMarkup(self:template(file['hydro/eqn/bssnok-fd.cl']))
end

function BSSNOKFiniteDifferenceEquationBase:initCodeModule_calcDT()
	self.solver.modules:add{
		name = 'calcDT',
		depends = table{
			'eqn.common',
			'coord_sqrt_g_ll##',
			'SETBOUNDS',
		}:append(
			self.cflMethod == '2008 Alcubierre' and { 
				'calc_gamma_uu',
			} or nil
		),
		code = self:template[[
kernel void calcDT(
	constant <?=solver.solver_t?>* solver,
	global real* dtBuf,
	const global <?=eqn.cons_t?>* UBuf,
	const global <?=solver.coord.cell_t?>* cellBuf
) {
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		dtBuf[index] = INFINITY;
		return;
	}
	real3 x = cellBuf[index].pos;
	const global <?=eqn.cons_t?>* U = UBuf + index;

<? if eqn.cflMethod == '2008 Alcubierre' then
?>	sym3 gamma_uu = calc_gamma_uu(U, x);
<? end 
?>
	real dt = INFINITY;
	<? for side=0,solver.dim-1 do ?>{
<? 
if eqn.cflMethod == '2013 Baumgarte et al, eqn 32' then
	-- TODO if the grid is static then this only needs to be done once
	if side == 0 then 
?>		dt = (real)min(dt, solver->grid_dx.x);
<?	elseif side == 1 then 
?>		dt = (real)min(dt, .5 * solver->grid_dx.x * solver->grid_dx.y);
<? 	elseif side == 2 then 
?>		dt = (real)min(dt, .5 * solver->grid_dx.x * sin(.5 * solver->grid_dx.y) * solver->grid_dx.z);
<? 	end 
else
	if eqn.cflMethod == '2008 Alcubierre' then 
?>		//this is asserting alpha and W >0, which they should be
		real absLambdaMax = U->alpha * sqrt(gamma_uu.<?=sym(side+1,side+1)?>);
		real dx = solver->grid_dx.s<?=side?>;
		dt = (real)min(dt, dx / absLambdaMax);
<? 	elseif eqn.cflMethod == '2017 Ruchlin et al, eqn 53' then 
?>		// for wavespeeds alpha sqrt(gammaBar^ii)
		// and if we assume alpha > 1
		// and gamma_ii ~ 1/gamma^ii 
		// and gammaBar_ij ~ gammaHat_ij 
		// then the typical CFL equation: dt <= dx / lambda, lambda = alpha sqrt(gammaBar^ii)
		// turns into the SENR code: dt <= sqrt(gammaHat_ii) * dx
		real sqrt_gammaHat_ii = coord_sqrt_g_ll<?=side..side?>(x);
		real ds = sqrt_gammaHat_ii  * solver->grid_dx.s<?=side?>;
		dt = (real)min(dt, ds);
<? 	else
		error("unknown cflMethod "..tostring(eqn.cflMethod))
	end
end
?>	}<? end ?>
	dtBuf[index] = dt;
}
]],
	}
end

function BSSNOKFiniteDifferenceEquationBase:getModuleDependsApplyInitCond()
	return {
		'setFlatSpace',
	}
end

return BSSNOKFiniteDifferenceEquationBase 

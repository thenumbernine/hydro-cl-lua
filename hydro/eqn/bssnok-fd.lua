local class = require 'ext.class'
local table = require 'ext.table'
local EinsteinEqn = require 'hydro.eqn.einstein'
local makePartials = require 'hydro.eqn.makepartial'

local common = require 'hydro.common'
local xNames = common.xNames


local BSSNOKFiniteDifferenceEquationBase = class(EinsteinEqn)

-- seems all the hyperbolic formalisms listed in Alcubierre's book use alpha sqrt(gamma^ii) for the speed-of-light wavespeed
-- however the 2017 Ruchlin paper says to use gamma_ij
--BSSNOKFiniteDifferenceEquationBase.cflMethod = '2008 Alcubierre'
--BSSNOKFiniteDifferenceEquationBase.cflMethod = '2013 Baumgarte et al, eqn 32'
BSSNOKFiniteDifferenceEquationBase.cflMethod = '2017 Ruchlin et al, eqn 53'

-- never need initDerivs for the finite-difference solvers
-- since they don't use any first-derivative state variables like the hyperbolic conservation law solvers do
BSSNOKFiniteDifferenceEquationBase.needsInitDerivs = false

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
		[8] = {[0] = -25, 48, -36, 16, -3, denom=12},
		senr_8 = { 	-- SENR's coeffs
			[-1] = -3,
			[0] = -10,
			[1] = 18,
			[2] = -6,
			[3] = 1,
			denom = 12,
		},
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
	self:initCodeModules_calc_gamma()
end

function BSSNOKFiniteDifferenceEquationBase:initCodeModules_calc_gamma()
	local solver = self.solver
	
	-- gammaHat_ij and co
	
	solver.modules:add{
		name = 'calc_gammaHat_ll',
		depends = {'coord_g_ll'},
		headercode = [[
#define calc_gammaHat_ll	coord_g_ll
]],
	}

	solver.modules:add{
		name = 'calc_gammaHat_uu',
		depends = {'coord_g_uu'},
		headercode = [[
#define calc_gammaHat_uu 	coord_g_uu
]],
	}

	solver.modules:add{
		name = 'calc_det_gammaHat',
		depends = {'coord_det_g'},
		code = [[
#define calc_det_gammaHat 	coord_det_g
]],
	}

	-- gammaBar_IJ and co

	solver.modules:add{
		name = 'calc_gammaHat_LL',
		headercode = [[
#define calc_gammaHat_LL(x) (sym3_ident)
]],
	}

	solver.modules:add{
		name = 'calc_gammaHat_UU',
		headercode = [[
#define calc_gammaHat_UU(x) (sym3_ident)
]],
	}

	solver.modules:add{
		name = 'calc_gammaBar_LL',
		depends = {
			'cons_t',
			'calc_gammaHat_LL',
		},
		code = self:template[[
sym3 calc_gammaBar_LL(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaHat_LL = calc_gammaHat_LL(x);
	sym3 gammaBar_LL = sym3_add(gammaHat_LL, U->epsilon_LL);
	return gammaBar_LL;
}
]],
	}

	solver.modules:add{
		name = 'calc_det_gammaBarLL',
		code = [[
/*
det(epsilon_IJ + gammaHat_IJ) 
= det(epsilon_IJ + delta_IJ) 
= det(e^i_I e^j_J (gammaHat_ij + epsilon_ij))
= det(e^i_I) det(e^j_J) det(gammaHat_ij + epsilon_ij)
= det(gammaHat^ij) det(gammaBar_ij)
= det(gammaBar_ij) / det(gammaHat_ij)
= 1
TODO detg ... unless we want to change the constraint
*/
#if 0	//use the value
real calc_det_gammaBarLL(global const <?=eqn.cons_t?>* U, ral3 x) {
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = sym3_det(gammaBar_LL);
	return det_gammaBarLL;
}
#else	//use the constraint
#define calc_det_gammaBarLL(x) 1.
#endif
]],
	}

	solver.modules:add{
		name = 'calc_gammaBar_UU',
		depends = {
			'cons_t',
			'calc_gammaBar_LL',
			'calc_det_gammaBarLL',
		},
		code = self:template[[
sym3 calc_gammaBar_UU(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	return gammaBar_UU;
}
]],
	}
	
	-- gammaBar_ij and co

	solver.modules:add{
		name = 'calc_gammaBar_ll',
		depends = {
			'cons_t',
			'calc_gammaHat_ll',
			'rescaleFromCoord/rescaleToCoord',
		},
		code = self:template[[
//gammaBar_ll.ij := gammaBar_ij = gammaHat_ij + epsilon_ij = gammaHat_ij + epsilon_IJ e_i^I e_j^J
sym3 calc_gammaBar_ll(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaHat_ll = calc_gammaHat_ll(x);
	sym3 epsilon_ll = sym3_rescaleToCoord_LL(U->epsilon_LL, x);
	sym3 gammaBar_ll = sym3_add(gammaHat_ll, epsilon_ll);
	return gammaBar_ll;
}
]],
	}

	solver.modules:add{
		name = 'calc_det_gammaBar',
		depends = {
			'calc_det_gammaHat',
		},
		code = self:template[[
//det(gammaBar_ij) = det(gammaHat_ij + epsilon_ij)
//...except sometimes, according to 2012 Baumgarte et al, last paragraph of II B
real calc_det_gammaBar(real3 x) {
	//TODO detg ...
	real det_gammaHat = calc_det_gammaHat(x);
	real detg = 1.;
	real det_gammaBar = det_gammaHat * detg;
	return det_gammaBar;
}
]],
	}

	solver.modules:add{
		name = 'calc_exp_neg4phi',
		headercode = [[
#define calc_exp_neg4phi(U) ((U)->W * (U)->W)
]],
	}

	solver.modules:add{
		name = 'calc_gammaBar_uu',
		depends = {
			'cons_t',
			'calc_gammaBar_ll',
			'calc_det_gammaBar',
		},
		code = self:template[[
sym3 calc_gammaBar_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
	return gammaBar_uu;
}
]],
	}

	solver.modules:add{
		name = 'calc_gamma_ll',
		depends = {
			'cons_t',
			'calc_gammaBar_ll',
			'calc_exp_neg4phi',
		},
		code = self:template[[
sym3 calc_gamma_ll(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	return gamma_ll;
}
]],
	}

	solver.modules:add{
		name = 'calc_gamma_uu',
		depends = {
			'cons_t',
			'calc_gammaBar_ll',
			'calc_exp_neg4phi',
			'calc_det_gammaBar',
		},
		code = self:template[[
sym3 calc_gamma_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	real det_gamma = calc_det_gammaBar(x) * exp_4phi * exp_4phi * exp_4phi;
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma); 
	return gamma_uu;
}
]],
	}



end

function BSSNOKFiniteDifferenceEquationBase:initCodeModule_setFlatSpace()
	local solver = self.solver
	solver.modules:add{
		name = 'mystery_C_U',
		code = [[
//TODO 2017 Ruchlin eqn. 8, what is C^i?
#define mystery_C_U	real3_zero
]],
	}

	solver.modules:add{
		name = 'setFlatSpace',
		depends = {'mystery_C_U'},
		code = self:template[[
void setFlatSpace(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* U,
	real3 x
) {
	U->alpha = 1.;
	U->beta_U = real3_zero;
	U->epsilon_LL = sym3_zero;
	U->W = 1;
	U->K = 0;
	U->ABar_LL = sym3_zero;

	//LambdaBar^i = Delta^i + C^i = Delta^i_jk gammaBar^jk = (connBar^i_jk - connHat^i_jk) gammaBar^jk + C^i
	//but when space is flat we have connBar^i_jk = connHat^i_jk and therefore Delta^i_jk = 0, Delta^i = 0, and LambdaBar^i = 0
	U->LambdaBar_U = mystery_C_U;

<? if eqn.useShift == 'HyperbolicGammaDriver' then
?>	U->B_U = real3_zero;
<? end
?>

	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0;
	U->M_U = real3_zero;
}
]],
	}
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

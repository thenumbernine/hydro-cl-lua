--[[
Baumgarte & Shapiro "Numerical Relativity: Solving Einstein's Equations on the Computer" 2010
Alcubierre "Introduction to Numerical Relativity" 2008

then I'm applying 2017 Ruchlin changes...
*) separate gammaBar_ll - gammaHat_ll = epsilon_ll
*) coordinate-transform beta^i, epsilon_ij, ABar_ij, LambdaBar^i to eliminate singularities from the metric

tensors are denoted with suffixes _u _l etc for upper and lower
rescaled tensors are denoted _U _L etc

Then I'm double checking all against (and borrowing heavily from) Zach Etienne's SENR: https://math.wvu.edu/~zetienne/SENR/
--]]
local file = require 'ext.file'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local EinsteinEqn = require 'eqn.einstein'
local makestruct = require 'eqn.makestruct'
local applyCommon = require 'common'
local time, getTime = table.unpack(require 'util.time')

local makePartials = require 'eqn.makepartial'

local BSSNOKFiniteDifferenceEquation = class(EinsteinEqn)
BSSNOKFiniteDifferenceEquation.name = 'BSSNOK finite difference' 
BSSNOKFiniteDifferenceEquation.hasEigenCode = true
BSSNOKFiniteDifferenceEquation.hasCalcDTCode = true
BSSNOKFiniteDifferenceEquation.hasFluxFromConsCode = true
BSSNOKFiniteDifferenceEquation.useConstrainU = true
BSSNOKFiniteDifferenceEquation.useSourceTerm = true

-- not used with finite-difference schemes anyways
BSSNOKFiniteDifferenceEquation.weightFluxByGridVolume = false

--[[
args:
	useShift = 'none'
				'GammaDriver'
				'HyperbolicGammaDriver' (default)
--]]
function BSSNOKFiniteDifferenceEquation:init(args)
	-- options:
	-- needs to be defined up front
	-- otherwise rebuild intVars based on it ...
	self.useShift = args.useShift or 'HyperbolicGammaDriver'

	local intVars = table{
		{name='alpha', type='real'},			-- 1
		{name='beta_U', type='real3'},		 	-- 3: beta^i
		{name='epsilon_LL', type='sym3'},		-- 6: gammaBar_ij - gammaHat_ij, only 5 dof since det gammaBar_ij = 1
		{name='W', type='real'},				-- 1: W = exp(-2 phi) = (det gammaHat_ij / det gamma_ij)^(1/6)
		{name='K', type='real'},				-- 1: K = K^i_i
		{name='ABar_LL', type='sym3'},			-- 6: ABar_ij, only 5 dof since ABar^k_k = 0
		{name='LambdaBar_U', type='real3'},		-- 3: LambdaBar^i = C^i + Delta^i = C^i + gammaBar^jk (connBar^i_jk - connHat^i_jk)
												-- TODO what is C^i ?
	}
	if self.useShift == 'HyperbolicGammaDriver' then
		intVars:insert{name='B_U', type='real3'}
	end

	self.consVars = table()
	:append(intVars)
	:append{
		--hyperbolic variables:
		--real3 a;				//3: a_i
		--_3sym3 dBar;			//18: dBar_ijk, only 15 dof since dBar_ij^j = 0
		--real3 Phi;			//3: Phi_i

		--stress-energy variables:
		{name='rho', type='real'},				--1: n_a n_b T^ab
		{name='S_u', type='real3'},				--3: -gamma^ij n_a T_aj
		{name='S_ll', type='sym3'},				--6: gamma_i^c gamma_j^d T_cd

		--constraints:
		{name='H', type='real'},				--1
		{name='M_u', type='real3'},				--3
	}
	self.numIntStates = makestruct.countScalars(intVars)
	
	-- call construction / build structures	
	BSSNOKFiniteDifferenceEquation.super.init(self, args)
end

function BSSNOKFiniteDifferenceEquation:createInitState()
	BSSNOKFiniteDifferenceEquation.super.createInitState(self)
	self:addGuiVars{
		{name='constrain_det_gammaBar', value=true, compileTime=true},
		--{name='constrain_det_gammaBar', value=false, compileTime=true},

		--{name='constrain_tr_ABar', value=true, compileTime=true},
		{name='constrain_tr_ABar', value=false, compileTime=true},
		
		{name='calc_H_and_M', value=true, compileTime=true},
		{name='diffuseSigma', value=.01},
		
		--{name='alphaMin', value=1e-7},
		{name='alphaMin', value=0},

		{name='shift_eta', value=1},	--1, or 1 / (2 M), for total mass M
	}
end

local string = require 'ext.string'

local function deduceRescale(name, fieldType)
	local parts = string.split(name, '_')
	-- if there's no _ in the variable name then it's probably a scalar
	assert((#parts == 1) == (fieldType == 'real'), "variable "..name.." is of type "..fieldType..", but I think it should be a real")
	if #parts == 1 then return false, parts end
	return parts[#parts], parts
end

--[[
copied and modified from eqn/makepartial.lua
I'm not sure if this will be my winning bssnok scheme
until then I'll keep the changes local

rescale = nil to not rescale, or the tensor suffix to use rescaling
	real3: suffix can be u or l
	sym3: suffix can be uu or ll
	_3sym3: suffix can be uuu or lll
	(see coord/coord.lua for more on rescaling)
	Notice, rescale functions need the coordinate.  I'm assuming 'x' for now.
	Right now I'm just rescaling before calculating the partial, but not after, so the partial index is in coordinate form.
--]]

local derivCoeffs = makePartials.derivCoeffs
local clnumber = require 'cl.obj.number'

local function makePartial(order, solver, field, fieldType, nameOverride, rescale)
	local suffix = 'l'
	if not field:find'_' then suffix = '_' .. suffix end

	local function add(x,y) return fieldType..'_add('..x..', '..y..')' end
	local function sub(x,y) return fieldType..'_sub('..x..', '..y..')' end
	local function real_mul(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
	local zero = fieldType..'_zero'
	local name = nameOverride or ('partial_'..field..suffix)

	local parts
	if rescale == true then 
		rescale, parts = deduceRescale(field, fieldType) 
		-- lower the last _... of the name
		if rescale and not nameOverride then
			local a, b = name:match'^(.*)_(.-)$'
			if a and b then
				b = b:lower()
				name = a..'_'..b
			end
		end
	end
	
	local d1coeffs = assert(derivCoeffs[1][order], "couldn't find 1st derivative coefficients of order "..order)
	local lines = table{'\t'..fieldType..' '..name..'[3];\n'}
	for i,xi in ipairs(xNames) do
		local namei = name..'['..(i-1)..']'
		local expr = zero
		if i <= solver.dim then
			for j,coeff in ipairs(d1coeffs) do
				local UR = 'U['..j..' * solver->stepsize.'..xi..'].'..field
				local UL = 'U[-'..j..' * solver->stepsize.'..xi..'].'..field
				if rescale then
					UL = fieldType..'_rescaleToCoord_'..rescale..'('..UL..', x)'
					UR = fieldType..'_rescaleToCoord_'..rescale..'('..UR..', x)'
				end
				expr = add(expr, real_mul(sub(UR, UL), clnumber(coeff)))
			end
			expr = real_mul(expr, '1. / solver->grid_dx.'..xi)
		end
		lines:insert('\t'..namei..' = '..expr..';')
	end
	return lines:concat'\n'
end

local function makePartial2(order, solver, field, fieldType, nameOverride, rescale)
	local suffix = 'll'
	if not field:find'_' then suffix = '_' .. suffix end
	
	local function add(x,y) return fieldType..'_add('..x..', '..y..')' end
	local function sub(x,y) return fieldType..'_sub('..x..', '..y..')' end
	local function real_mul(x,y) return fieldType..'_real_mul('..x..', '..y..')' end
	local zero = fieldType..'_zero'
	local name = nameOverride or ('partial2_'..field..suffix)

	local parts
	if rescale == true then 
		rescale, parts = deduceRescale(field, fieldType) 
		-- lower the last _... of the name
		if rescale and not nameOverride then
			local a, b = name:match'^(.*)_(.-)$'
			if a and b then
				b = b:lower()
				name = a..'_'..b
			end	
		end
	end

	local d1coeffs = assert(derivCoeffs[1][order], "couldn't find 1st derivative coefficients of order "..order)
	local d2coeffs = assert(derivCoeffs[2][order], "couldn't find 2nd derivative coefficients of order "..order)
	local lines = table()
	lines:insert('\t'..fieldType..' '..name..'[6];')
	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi, xj = xNames[i], xNames[j]
		local nameij = name..'['..(ij-1)..']'
		if i > solver.dim or j > solver.dim then
			lines:insert('\t'..nameij..' = '..zero..';')
		elseif i == j then
			local expr = real_mul('U->'..field, d2coeffs[0])
			for k,coeff in ipairs(d2coeffs) do
				local UR = 'U['..k..' * solver->stepsize.s'..(i-1)..'].'..field
				local UL = 'U[-'..k..' * solver->stepsize.s'..(i-1)..'].'..field
				if rescale then
					UL = fieldType..'_rescaleToCoord_'..rescale..'('..UL..', x)'
					UR = fieldType..'_rescaleToCoord_'..rescale..'('..UR..', x)'
				end
				expr = add(expr, real_mul(add(UR, UL), clnumber(coeff)))
			end
			expr = real_mul(expr, '1. / (solver->grid_dx.'..xi..' * solver->grid_dx.'..xi..')')
			lines:insert('\t'..nameij..' = '..expr..';')
		else
			local expr = zero
			for k,coeff_k in ipairs(d1coeffs) do
				for l,coeff_l in ipairs(d1coeffs) do
					local URR = 'U['..k..' * solver->stepsize.'..xi..' + '..l..' * solver->stepsize.'..xj..'].'..field
					local ULL = 'U[-'..k..' * solver->stepsize.'..xi..' - '..l..' * solver->stepsize.'..xj..'].'..field
					local ULR = 'U[-'..k..' * solver->stepsize.'..xi..' + '..l..' * solver->stepsize.'..xi..'].'..field
					local URL = 'U['..k..' * solver->stepsize.'..xi..' - '..l..' * solver->stepsize.'..xi..'].'..field
					if rescale then
						ULL = fieldType..'_rescaleToCoord_'..rescale..'('..ULL..', x)'
						ULR = fieldType..'_rescaleToCoord_'..rescale..'('..ULR..', x)'
						URL = fieldType..'_rescaleToCoord_'..rescale..'('..URL..', x)'
						URR = fieldType..'_rescaleToCoord_'..rescale..'('..URR..', x)'
					end
					expr = add(expr, real_mul(sub(add(URR, ULL), add(ULR, URL)), clnumber(coeff_k * coeff_l)))
				end
			end
			expr = real_mul(expr, '1. / (solver->grid_dx.'..xi..' * solver->grid_dx.'..xj..')')
			lines:insert('\t'..nameij..' = '..expr..';')
		end
	end
	return lines:concat'\n'
end

function BSSNOKFiniteDifferenceEquation:makePartial(field, fieldType, nameOverride, rescale)
	local derivOrder = 2 * self.solver.numGhost
	if fieldType == nil then
		local _, var = self.consVars:find(nil, function(v) return v.name == field end)
		assert(var)
		fieldType = var.type
	end
	if rescale == nil then rescale = true end
	return makePartial(derivOrder, self.solver, field, fieldType, nameOverride, rescale)
end

function BSSNOKFiniteDifferenceEquation:makePartial2(field, fieldType, nameOverride, rescale)
	local derivOrder = 2 * self.solver.numGhost
	if fieldType == nil then
		local _, var = self.consVars:find(nil, function(v) return v.name == field end)
		assert(var)
		fieldType = var.type
	end
	if rescale == nil then rescale = true end
	return makePartial2(derivOrder, self.solver, field, fieldType, nameOverride, rescale)
end

function BSSNOKFiniteDifferenceEquation:compile(expr)
	local symmath = require 'symmath'
	local baseCoords = self.solver.coord.baseCoords
	local var = symmath.var
	return symmath.export.C(expr
		:replace(baseCoords[1], var'x.x')
		:replace(baseCoords[2], var'x.y')
		:replace(baseCoords[3], var'x.z'))
end

function BSSNOKFiniteDifferenceEquation:getEnv()
	return applyCommon{
		eqn = self,
		solver = self.solver,
	}
end

function BSSNOKFiniteDifferenceEquation:getCommonFuncCode()
	return template([[

<? local dim = solver.dim ?>

#if 1	
//rescaling, used for bssn finite-difference, but I am tempted to try it with other coordinate systems with singularities
//TODO for the initial conditions do this symbolically instead of numerically

//apply this to lower indexes to convert from coordinate metric to better metric
//apply this to upper indexes to convert from better metric to coordinate metric
real3 real3_rescaleFromCoord_l(real3 v, real3 x) {
	return (real3){
		.x = v.x / coord_dx0(x),
		.y = v.y / coord_dx1(x),
		.z = v.z / coord_dx2(x),
	};
}
#define real3_rescaleToCoord_U real3_rescaleFromCoord_l

//convert coord upper to better
//convert better lower to coord
real3 real3_rescaleToCoord_L(real3 v, real3 x) {
	return (real3){
		.x = v.x * coord_dx0(x),
		.y = v.y * coord_dx1(x),
		.z = v.z * coord_dx2(x),
	};
}
#define real3_rescaleFromCoord_u real3_rescaleToCoord_L

sym3 sym3_rescaleFromCoord_ll(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
<? end
?>	};
}
#define sym3_rescaleToCoord_UU sym3_rescaleFromCoord_ll

sym3 sym3_rescaleToCoord_LL(sym3 a, real3 x) {
	return (sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = a.<?=xij?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x)),
<? end
?>	};
}
#define sym3_rescaleFromCoord_uu sym3_rescaleToCoord_LL

_3sym3 _3sym3_rescaleFromCoord_lll(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define _3sym3_rescaleToCoord_UUU _3sym3_rescaleFromCoord_lll

_3sym3 _3sym3_rescaleToCoord_LLL(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define _3sym3_rescaleFromCoord_uuu _3sym3_rescaleToCoord_LLL


_3sym3 _3sym3_rescaleFromCoord_ull(_3sym3 a, real3 x) {
	return (_3sym3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = (sym3){
<?	for jk,xjk in ipairs(symNames) do
	local j,k = from6to3x3(jk)
?>			a.<?=xi?>.<?=xjk?> * coord_dx<?=i-1?>(x) / (coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x)),
<?	end
?>		},
<? end
?>	};
}


sym3sym3 sym3sym3_rescaleFromCoord_lll(sym3sym3 a, real3 x) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = (sym3){
<?	for kl,xkl in ipairs(symNames) do
	local k,l = from6to3x3(kl)
?>			.<?=xkl?> = a.<?=xij?>.<?=xkl?> / (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x) * coord_dx<?=l-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define sym3sym3_rescaleToCoord_UUUU sym3sym3_rescaleFromCoord_llll

sym3sym3 sym3sym3_rescaleToCoord_LLLL(sym3sym3 a, real3 x) {
	return (sym3sym3){
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		.<?=xij?> = (sym3){
<?	for kl,xkl in ipairs(symNames) do
	local k,l = from6to3x3(kl)
?>			.<?=xkl?> = a.<?=xij?>.<?=xkl?> * (coord_dx<?=i-1?>(x) * coord_dx<?=j-1?>(x) * coord_dx<?=k-1?>(x) * coord_dx<?=l-1?>(x)),
<?	end
?>		},
<? end
?>	};
}
#define sym3sym3_rescaleFromCoord_uuuu sym3sym3_rescaleToCoord_LLLL

#else	//debugging -- turning it off

#define real3_rescaleFromCoord_l(a,x) a
#define real3_rescaleToCoord_U(a,x) a
#define real3_rescaleToCoord_L(a,x) a
#define real3_rescaleFromCoord_u(a,x) a
#define sym3_rescaleFromCoord_ll(a,x) a
#define sym3_rescaleToCoord_UU(a,x) a
#define sym3_rescaleToCoord_LL(a,x) a
#define sym3_rescaleFromCoord_uu(a,x) a
#define _3sym3_rescaleFromCoord_lll(a,x) a
#define _3sym3_rescaleToCoord_UUU(a,x) a
#define _3sym3_rescaleToCoord_LLL(a,x) a
#define _3sym3_rescaleFromCoord_uuu(a,x) a
#define sym3sym3_rescaleFromCoord_lll(a,x) a
#define sym3sym3_rescaleToCoord_UUUU(a,x) a
#define sym3sym3_rescaleToCoord_LLLL(a,x) a
#define sym3sym3_rescaleFromCoord_uuuu (a,x) a

#endif


//TODO 2017 Ruchlin eqn. 8, what is C^i?
#define mystery_C_U	real3_zero


	// gammaHat_ij and co


#define calc_gammaHat_ll	coord_g_ll
#define calc_det_gammaHat 	coord_det_g
#define calc_gammaHat_uu 	coord_g_uu
#define calc_connHat_ull	coord_conn_ull
#define calc_connHat_lll	coord_conn_lll

//partial_gammaHat_lll.k.ij := gammaHat_ij,k
_3sym3 calc_partial_gammaHat_lll(real3 x) {
	_3sym3 partial_gammaHat_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
?>	partial_gammaHat_lll.<?=xk?>.<?=xij?> = <?=eqn:compile(solver.coord.dg[k][i][j])?>;
<?	end
end
?>	return partial_gammaHat_lll;
}

//partial2_gammaHat_llll.kl.ij := gammaHat_ij,kl
sym3sym3 calc_partial2_gammaHat_llll(real3 x) {
	sym3sym3 partial2_gammaHat_llll;
<?
local symmath = require 'symmath'
local Tensor = symmath.Tensor
local d2g = Tensor'_ijkl'
d2g['_klij'] = solver.coord.dg'_kij,l'()
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	for kl, xkl in ipairs(symNames) do
		local k,l = from6to3x3(kl)
?>	partial2_gammaHat_llll.<?=xij?>.<?=xkl?> = <?=eqn:compile(d2g[k][l][i][j])?>;
<?	end
end
?>	return partial2_gammaHat_llll;
}

//connHat^i_jk,l := partial_connHat_ulll[l].i.jk
void calc_partial_connHat_ulll(_3sym3 partial_connHat_ulll[3], real3 x) {
<? 
local dGamma = solver.coord.Gamma_ull'^i_jk,l'()
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
		for l,xl in ipairs(xNames) do
?>	partial_connHat_ulll[<?=l-1?>].<?=xi?>.<?=xjk?> = <?=eqn:compile(dGamma[i][j][k][l])?>;
<?		end
	end
end
?>}

real3 calc_partial_det_gammaHat_l(real3 x) {
	real3 partial_det_gammaHat_l;
<? 
local partial_det_gammaHat_l = Tensor('_i', function(i)
	return solver.coord.det_g:diff(solver.coord.baseCoords[i])()
end)
for i,xi in ipairs(xNames) do
?>	partial_det_gammaHat_l.<?=xi?> = <?=eqn:compile(partial_det_gammaHat_l[i])?>;
<? end
?>	return partial_det_gammaHat_l;
}

sym3 calc_partial2_det_gammaHat_ll(real3 x) {
	sym3 partial2_det_gammaHat_ll;
<?
local partial2_det_gammaHat_ll = partial_det_gammaHat_l'_i,j'():factorDivision()
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>	partial2_det_gammaHat_ll.<?=xij?> = <?=eqn:compile(partial2_det_gammaHat_ll[i][j])?>;
<?
end
?>	return partial2_det_gammaHat_ll;
}


	// gammaBar_IJ and co


static inline sym3 calc_gammaHat_LL(real3 x) {
	return sym3_ident;
}

static inline sym3 calc_gammaHat_UU(real3 x) {
	return sym3_ident;
}

sym3 calc_gammaBar_LL(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaHat_LL = calc_gammaHat_LL(x);
	sym3 gammaBar_LL = sym3_add(gammaHat_LL, U->epsilon_LL);
	return gammaBar_LL;
}

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
#define calc_det_gammaBarLL(x) 1

sym3 calc_gammaBar_UU(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	return gammaBar_UU;
}


	// gammaBar_ij and co


//gammaBar_ll.ij := gammaBar_ij = gammaHat_ij + epsilon_ij = gammaHat_ij + epsilon_IJ e_i^I e_j^J
sym3 calc_gammaBar_ll(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaHat_ll = calc_gammaHat_ll(x);
	sym3 epsilon_ll = sym3_rescaleToCoord_LL(U->epsilon_LL, x);
	sym3 gammaBar_ll = sym3_add(gammaHat_ll, epsilon_ll);
	return gammaBar_ll;
}

//det(gammaBar_ij) = det(gammaHat_ij + epsilon_ij)
//...except sometimes, according to 2012 Baumgarte et al, last paragraph of II B
real calc_det_gammaBar(real3 x) {
	//TODO detg ...
	real det_gammaHat = calc_det_gammaHat(x);
	real detg = 1.;
	real det_gammaBar = det_gammaHat * detg;
	return det_gammaBar;
}

#define calc_exp_neg4phi(U) ((U)->W * (U)->W)

sym3 calc_gammaBar_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
	return gammaBar_uu;
}

sym3 calc_gamma_uu(global const <?=eqn.cons_t?>* U, real3 x) {
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	real det_gamma = calc_det_gammaBar(x) * exp_4phi * exp_4phi * exp_4phi;
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma); 
	return gamma_uu;
}

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
	U->M_u = real3_zero;
}

]], self:getEnv())
end

function BSSNOKFiniteDifferenceEquation:getCode_connBar_ull()
 	return template([[
	_3sym3 connBar_ull;
	{
		//connBar_lll.i.jk := connBar_ijk = 1/2 (gammaBar_ij,k + gammaBar_ik,j - gammaBar_jk,i)
		_3sym3 connBar_lll;
<? 
for i,xi in ipairs(xNames) do
	for jk,xjk in ipairs(symNames) do
		local j,k = from6to3x3(jk)
		local xj,xk = xNames[j],xNames[k]
?>		connBar_lll.<?=xi?>.<?=xjk?> = .5 * (0.
			+ partial_gammaBar_lll.<?=xk?>.<?=sym(i,j)?>
			+ partial_gammaBar_lll.<?=xj?>.<?=sym(i,k)?>
			- partial_gammaBar_lll.<?=xi?>.<?=sym(j,k)?>
		);
<?	end
end
?>		//connBar_ull[i].jk := connBar^i_jk = gammaBar^il connBar_ljk
		connBar_ull = sym3_3sym3_mul(gammaBar_uu, connBar_lll);
	}
]], self:getEnv())
end

function BSSNOKFiniteDifferenceEquation:getCode_RBar_LL()
	return template([[

	//partial2_gammaBar_llll.kl.ij := gammaBar_ij,kl
	sym3sym3 partial2_gammaBar_llll;
<? 
for ij,xij in ipairs(symNames) do
	for kl,xkl in ipairs(symNames) do
?>	partial2_gammaBar_llll.<?=xkl?>.<?=xij?> = partial2_epsilon_llll[<?=kl-1?>].<?=xij?> + partial2_gammaHat_llll.<?=xkl?>.<?=xij?>;
<?	end
end 
?>

	//DHat_gammaBar_lll.k.ij = DHat_k gammaBar_ij 
	// = gammaBar_ij,k - connHat^l_ki gammaBar_lj - connHat^l_kj gammaBar_il
	_3sym3 DHat_gammaBar_lll;
<?
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
		local i,j = from6to3x3(ij)
		local xi,xj = xNames[i], xNames[j]
?>	DHat_gammaBar_lll.<?=xk?>.<?=xij?> = partial_gammaBar_lll.<?=xk?>.<?=xij?>
<?		for l,xl in ipairs(xNames) do
?>		- connHat_ull.<?=xl?>.<?=sym(k,i)?> * gammaBar_ll.<?=sym(l,j)?>
		- connHat_ull.<?=xl?>.<?=sym(k,j)?> * gammaBar_ll.<?=sym(l,i)?>
<?		end
?>	;
<?	end
end
?>

	/*
	partial_DHat_gammaBar_llll[l].k.ij := partial_l DHat_k gammaBar_ij
		= (gammaBar_ij,k - connHat^m_ki gammaBar_mj - connHat^m_kj gammaBar_mi)_,l
		= gammaBar_ij,kl 
			- connHat^m_ki,l gammaBar_mj 
			- connHat^m_kj,l gammaBar_mi
			- connHat^m_ki gammaBar_mj,l
			- connHat^m_kj gammaBar_mi,l
	*/
	_3sym3 partial_DHat_gammaBar_llll[3];
<? 
for k,xk in ipairs(xNames) do
	for l,xl in ipairs(xNames) do
		for ij,xij in ipairs(symNames) do
			local i,j = from6to3x3(ij)
			local xi,xj = xNames[i], xNames[j]
?>	partial_DHat_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?> = 0.
		+ partial2_gammaBar_llll.<?=sym(k,l)?>.<?=xij?>
<?			for m,xm in ipairs(xNames) do
?>
		- partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(k,i)?> * gammaBar_ll.<?=sym(m,j)?>		//diverging: RBar_ij = diag(0, -.5, -.75)
		- partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(k,j)?> * gammaBar_ll.<?=sym(m,i)?>		//diverging: RBar_ij = diag(0, -.5, -.75)
		- connHat_ull.<?=xm?>.<?=sym(k,i)?> * partial_gammaBar_lll.<?=xl?>.<?=sym(m,j)?>		//diverging: RBar_ij = diag(0, 1, 1)
		- connHat_ull.<?=xm?>.<?=sym(k,j)?> * partial_gammaBar_lll.<?=xl?>.<?=sym(m,i)?>		//diverging: RBar_ij = diag(0, 1, 1)
<?			end
?>	;
<?		end
	end
end
?>

	/*
	DHat2_gammaBar_llll[l].k.ij = DHat_l DHat_k gammaBar_ij
		= partial_l DHat_k gammaBar_ij
			- connHat^m_lk DHat_m gammaBar_ij
			- connHat^m_li DHat_k gammaBar_mj
			- connHat^m_lj DHat_k gammaBar_im
	*/
	_3sym3 DHat2_gammaBar_llll[3];
<?
for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i], xNames[j]
	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>	DHat2_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?> = 0.
		+ partial_DHat_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?>
<?			for m,xm in ipairs(xNames) do
?>		- connHat_ull.<?=xm?>.<?=sym(l,k)?> * DHat_gammaBar_lll.<?=xm?>.<?=sym(i,j)?>
		- connHat_ull.<?=xm?>.<?=sym(l,i)?> * DHat_gammaBar_lll.<?=xk?>.<?=sym(m,j)?>
		- connHat_ull.<?=xm?>.<?=sym(l,j)?> * DHat_gammaBar_lll.<?=xk?>.<?=sym(i,m)?>
<?			end
?>	;
<?		end
	end
end
?>

	sym3 trBar_partial2_gammaBar_ll;
<? for ij,xij in ipairs(symNames) do
?>	trBar_partial2_gammaBar_ll.<?=xij?> = 0.
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * partial2_gammaBar_llll.<?=sym(k,l)?>.<?=xij?>
<?		end
	end
?>	;
<? end
?>

	//trBar_DHat2_gammaBar_ll.ij := gammaBar^kl DHat_k DHat_l gammaBar_ij
	sym3 trBar_DHat2_gammaBar_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	trBar_DHat2_gammaBar_ll.<?=xij?> = 0.
		+ trBar_partial2_gammaBar_ll.<?=xij?>
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * DHat2_gammaBar_llll[<?=l-1?>].<?=xk?>.<?=xij?>
<?		end
	end
?>	;
<? end
?>

	//derivative is the last index, unlike the partial_*'s
	//DHat_LambdaBar_ul.i.j := DHat_j LambdaBar^i = LambdaBar^i_,j + connHat^i_jk LambdaBar^k
	real3x3 DHat_LambdaBar_UL;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	DHat_LambdaBar_UL.<?=xi?>.<?=xj?> = (0.
		+ partial_LambdaBar_ul[<?=j-1?>].<?=xi?>
			 * coord_dx<?=i-1?>(x) / coord_dx<?=j-1?>(x))
<?		for k,xk in ipairs(xNames) do
?>		+ connHat_ULL.<?=xi?>.<?=sym(j,k)?> * U->LambdaBar_U.<?=xk?>
<?		end
?>	;
<?	end
end
?>

	/*
	2017 Ruchlin eqn 12
	RBar_ij = 
		-1/2 gammaBar^kl DHat_k DHat_l gammaBar_ij
		+ 1/2 gammaBar_ki DHat_j LambdaBar^k
		+ 1/2 gammaBar_kj DHat_i LambdaBar^k
		+ 1/2 Delta^k Delta_ijk
		+ 1/2 Delta^k Delta_jik
		+ gammaBar^kl (
			Delta^m_ki Delta_jml
			+ Delta^m_kj Delta_iml
			+ Delta^m_ik Delta_mjl
		)
	
	RBar_ij = 
		-1/2 gammaBar^kl DHat_k DHat_l gammaBar_ij
		+ 1/2 gammaBar_ki DHat_j LambdaBar^k
		+ 1/2 gammaBar_kj DHat_i LambdaBar^k
		+ 1/2 Delta^k Delta_ikj
		+ 1/2 Delta^k Delta_jki
		+ Delta^m_ki Delta_jm^k
		+ Delta^m_kj Delta_im^k
		+ Delta^m_ik Delta_mj^k
	*/
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	RBar_LL.<?=xij?> = 0.
			- .5 * trBar_DHat2_gammaBar_ll.<?=xij?>	
<?	for k,xk in ipairs(xNames) do
?>			+ .5 * gammaBar_LL.<?=sym(i,k)?> * DHat_LambdaBar_UL.<?=xk?>.<?=xj?> 
			+ .5 * gammaBar_LL.<?=sym(j,k)?> * DHat_LambdaBar_UL.<?=xk?>.<?=xi?>
			+ .5 * Delta_U.<?=xk?> * (
				Delta_LLL.<?=xi?>.<?=sym(k,j)?> 
				+ Delta_LLL.<?=xj?>.<?=sym(k,i)?>
			)
<?		for l,xl in ipairs(xNames) do
			for m,xm in ipairs(xNames) do
?>			+ gammaBar_UU.<?=sym(k,l)?> * (0.
				+ Delta_ULL.<?=xm?>.<?=sym(k,i)?> * Delta_LLL.<?=xj?>.<?=sym(m,l)?>
				+ Delta_ULL.<?=xm?>.<?=sym(k,j)?> * Delta_LLL.<?=xi?>.<?=sym(m,l)?>
				+ Delta_ULL.<?=xm?>.<?=sym(i,k)?> * Delta_LLL.<?=xm?>.<?=sym(j,l)?>
			)
<?			end
		end
	end
?>
	;
<? end ?>
]], self:getEnv())
end

function BSSNOKFiniteDifferenceEquation:getInitStateCode()
	return template([=[
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	real3 xc = coordMap(x);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	sym3 gammaHat_ll = calc_gammaHat_ll(x);

	real alpha = 1.;
	real3 beta_u = real3_zero;
	sym3 gamma_ll = gammaHat_ll;	
	sym3 K_ll = sym3_zero;
	real rho = 0.;

<? if eqn.initState.name == 'Minkowski' then ?>
	
	setFlatSpace(solver, U, x);

<? else -- not Minkowski ?>

	<?=code?>

	U->alpha = alpha;
	U->beta_U = real3_rescaleFromCoord_u(beta_u, x);

	real det_gamma = sym3_det(gamma_ll);
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma);
	
	//det(gammaBar_ij) == det(gammaHat_ij)
	real det_gammaBar = calc_det_gammaBar(x); 

	//gammaBar_ij = e^(-4phi) gamma_ij
	//real exp_neg4phi = exp(-4 * U->phi);
	real exp_neg4phi = cbrt(det_gammaBar / det_gamma);

	//W = exp(-2 phi)
	U->W = sqrt(exp_neg4phi);

	sym3 gammaBar_ll = sym3_real_mul(gamma_ll, exp_neg4phi);
	sym3 epsilon_ll = sym3_sub(gammaBar_ll, gammaHat_ll);
	U->epsilon_LL = sym3_rescaleFromCoord_ll(epsilon_ll, x);

	U->K = sym3_dot(K_ll, gamma_uu);
	sym3 A_ll = sym3_sub(K_ll, sym3_real_mul(gamma_ll, 1./3. * U->K));
	sym3 ABar_ll = sym3_real_mul(A_ll, exp_neg4phi);
	U->ABar_LL = sym3_rescaleFromCoord_ll(ABar_ll, x);

<? end -- Minkowski ?>

	U->rho = rho;
	U->S_u = real3_zero;
	U->S_ll = sym3_zero;
	U->H = 0.;
	U->M_u = real3_zero;
}

//after popularing gammaBar_ll, use its finite-difference derivative to initialize LambdaBar_u
//TODO do this symbolically.  That's what I originally did, but symbolic calculations were getting complex
// however, with spherical BSSN, you need to 
kernel void initDerivs(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;

#if 0	
	SETBOUNDS(0,0);
	if (OOB(numGhost,numGhost)) {
		setFlatSpace(solver, U, x);
		return;
	}
#endif
 
 <?=eqn:makePartial'epsilon_LL'?>

	_3sym3 partial_gammaHat_lll = calc_partial_gammaHat_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>
	
	_3sym3 connHat_lll = calc_connHat_lll(x);
	_3sym3 connHat_ull = calc_connHat_ull(x);
	sym3 gammaBar_uu = calc_gammaBar_uu(U, x);
	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);

<?=eqn:getCode_connBar_ull()?>
	
	//Delta^i_jk = connBar^i_jk - connHat^i_jk
	_3sym3 Delta_ULL = _3sym3_rescaleFromCoord_ull(_3sym3_sub(connBar_ull, connHat_ull), x);

	real3 LambdaBar_u;
	U->LambdaBar_U = _3sym3_sym3_dot23(Delta_ULL, gammaBar_UU);

}
]=], self:getEnv())
end

BSSNOKFiniteDifferenceEquation.solverCodeFile = 'eqn/bssnok-fd-num.cl'

function BSSNOKFiniteDifferenceEquation:getEigenTypeCode()
	return template([[
typedef struct { char unused; } <?=eqn.eigen_t?>;
]], {eqn=self})
end

BSSNOKFiniteDifferenceEquation.predefinedDisplayVars = {
--[=[
	'U alpha',
--	'U beta_U mag',
	'U beta_U x',
	'U beta_U y',
	'U beta_U z',
--	'U B_U mag',
	'U B_U x',
	'U B_U y',
	'U B_U z',
	--'U epsilon_LL norm',
	'U epsilon_LL xx',
	'U epsilon_LL xy',
	'U epsilon_LL xz',
	'U epsilon_LL yy',
	'U epsilon_LL yz',
	'U epsilon_LL zz',
	'U W',
	'U K',
--	'U ABar_LL tr weighted',
--] =]	
	'U ABar_LL xx',
	'U ABar_LL xy',
	'U ABar_LL xz',
	'U ABar_LL yy',
	'U ABar_LL yz',
	'U ABar_LL zz',
-- [=[	
--	'U LambdaBar_U mag',
	'U LambdaBar_U x',
	'U LambdaBar_U y',
	'U LambdaBar_U z',
	'U H',
--	'U M_u mag',
	'U M_u x',
	'U M_u y',
	'U M_u z',
	--'U det gammaBar - det gammaHat',
	--'U det gamma_ij based on phi',
	--'U volume',
	--'U f',
	--'U gamma_LL tr weighted',
--]=]	

-- [[ debugging derivatives
--[=[	
	'deriv alpha',
	'deriv beta_U x',
	'deriv beta_U y',
	'deriv beta_U z',
	'deriv B_U x',
	'deriv B_U y',
	'deriv B_U z',
	'deriv epsilon_LL xx',
	'deriv epsilon_LL xy',
	'deriv epsilon_LL xz',
	'deriv epsilon_LL yy',
	'deriv epsilon_LL yz',
	'deriv epsilon_LL zz',
	'deriv W',
	'deriv K',
--] =]	
	'deriv ABar_LL xx',
	'deriv ABar_LL xy',
	'deriv ABar_LL xz',
	'deriv ABar_LL yy',
	'deriv ABar_LL yz',
	'deriv ABar_LL zz',
-- [=[	
	'deriv LambdaBar_U x',
	'deriv LambdaBar_U y',
	'deriv LambdaBar_U z',
	'deriv H',
	'deriv M_u x',
	'deriv M_u y',
	'deriv M_u z',
--]=]	
--]]

	--'U tr_DBar2_phi',
	--'U DBar_phi_sq',
	--'U ABarSq tr weighted',

-- [[
	'U RBar_LL xx',
	'U RBar_LL xy',
	'U RBar_LL xz',
	'U RBar_LL yy',
	'U RBar_LL yz',
	'U RBar_LL zz',
--]]
--[[
	'U del gammaBar_ll sym xx',
	'U del gammaBar_ll sym xy',
	'U del gammaBar_ll sym xz',
	'U del gammaBar_ll sym yy',
	'U del gammaBar_ll sym yz',
	'U del gammaBar_ll sym zz',
--]]
}

function BSSNOKFiniteDifferenceEquation:getDisplayVars()	
	local vars = BSSNOKFiniteDifferenceEquation.super.getDisplayVars(self)

	vars:append{
		{
			name = 'gamma_ll',
			type = 'sym3',
			code = [[
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, exp_4phi);
	*value_sym3 = gamma_ll;
]], 
		},
		{name='gammaHat_ll', code=[[	*value_sym3 = calc_gammaHat_ll(x);]], type='sym3'},
		{name='gammaBar_ll', code=[[	*value_sym3 = calc_gammaBar_ll(U, x);]], type='sym3'},
		{name='gamma_uu', code=[[	*value_sym3 = calc_gamma_uu(U, x);]], type='sym3'},
		{name='gammaHat_uu', code=[[	*value_sym3 = calc_gammaHat_uu(x);]], type='sym3'},
		{name='gammaBar_uu', code=[[	*value_sym3 = calc_gammaBar_uu(U, x);]], type='sym3'},
		{name='K_ll', code=[[
	real exp_4phi = 1. / calc_exp_neg4phi(U);
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	*value_sym3 = sym3_real_mul(
		sym3_add(
			sym3_rescaleToCoord_LL(U->ABar_LL, x),
			sym3_real_mul(gammaBar_ll, U->K / 3.)
		), exp_4phi);
]], type='sym3'},

		{name='det gammaBar - det gammaHat', code=[[
	*value = sym3_det(calc_gammaBar_ll(U, x)) - calc_det_gammaBar(x);
]]},
		{name='det gamma based on phi', code=[[
	real exp_neg4phi = calc_exp_neg4phi(U);
	real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	real det_gamma = exp_12phi * calc_det_gammaHat(x);
	*value = det_gamma;
]]},
		
		{name='S', code='*value = sym3_dot(U->S_ll, calc_gamma_uu(U, x));'},
		
		{
			name='volume', 
			code=[[
	//|g| = exp(12 phi) |g_grid|
	real exp_neg4phi = calc_exp_neg4phi(U);
	real exp_12phi = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	real det_gamma = exp_12phi * calc_det_gammaHat(x);
	*value = U->alpha * det_gamma;
]],
		},
		{name='f', code='*value = calc_f(U->alpha);'},
		{name='df/dalpha', code='*value = calc_dalpha_f(U->alpha);'},
	
--[=[	
		{
			name = 'ABarSq',
			type = 'sym3',
			code = [[
	sym3 gammaBar_uu = calc_gammaBar_uu(U, x);
	sym3 gammaBar_UU = sym3_rescaleFromCoord_uu(gammaBar_uu, x);
	real3x3 ABar_UL = sym3_sym3_mul(gammaBar_UU, U->ABar_LL);
	sym3 ABarSq_LL = sym3_real3x3_to_sym3_mul(U->ABar_LL, ABar_UL);
	*value_sym3 = ABarSq_LL;
]],
		},
		
		{	-- gammaBar^ij DBar_i DBar_j phi
			name = 'tr_DBar2_phi',
			code = template([[
	
	<?=eqn:makePartial'epsilon_LL'?>

	_3sym3 partial_gammaHat_lll = calc_partial_gammaHat_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>
	
	
	sym3 gammaBar_uu = calc_gammaBar_uu(U, x);

<?=eqn:getCode_connBar_ull()?>

<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
	
	real3 partial_phi_l;
	sym3 partial2_phi_ll;
	{
		
		//partial_phi_l.i := phi_,i = -W_,i / (2 W) 
<? for i,xi in ipairs(xNames) do
?>		partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>

		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				-partial2_W_ll[<?=ij-1?>] 
				+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
			) / U->W;
<? end ?>
	}

	real tr_DBar2_phi = 0.
<? for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
		local ij = from3x3to6(i,j)
		local xij = symNames[ij]
?>		+ gammaBar_uu.<?=sym(i,j)?> * (
			partial2_phi_ll.<?=xij?>
<?		for k,xk in ipairs(xNames) do
?>			- connBar_ull.<?=xk?>.<?=xij?> * partial_phi_l.<?=xk?>
<?		end
?>		)
<?	end
end
?>	;

	*value = tr_DBar2_phi;
]], self:getEnv())
		},


		{
			name = 'DBar2_phi_ll',
			code = template([[
<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
	
	real3 partial_phi_l;
	sym3 partial2_phi_ll;
	{
		
		//partial_phi_l.i := phi_,i = -W_,i / (2 W) 
<? for i,xi in ipairs(xNames) do
?>		partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>

		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				-partial2_W_ll[<?=ij-1?>] 
				+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
			) / U->W;
<? end ?>
	}

<?=eqn:makePartial'epsilon_LL'?>
	_3sym3 partial_gammaHat_lll = calc_partial_gammaHat_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>

	sym3 gammaBar_uu = calc_gammaBar_uu(U, x);
<?=eqn:getCode_connBar_ull()?>
	
	//DBar2_phi_ll.ij := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	sym3 DBar2_phi_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_ll.<?=xij?> = partial2_phi_ll.<?=xij?> 
<?	for k,xk in ipairs(xNames) do	
?>		- connBar_ull.<?=xk?>.<?=xij?> * partial_phi_l.<?=xk?>
<?	end
?>	;
<? end
?>
	*value = sym3_dot(gammaBar_uu, DBar2_phi_ll);
]], self:getEnv()),
 		},
 
		{
			name = 'partial_phi_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial'W'?>
	real3 partial_phi_l;
<? for i,xi in ipairs(xNames) do
?>	partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>
	
	*value_real3 = partial_phi_l;
]], self:getEnv()),
		},

		{
			name = 'partial_alpha_l',
			type = 'real3',
			code = template([[
<?=eqn:makePartial'alpha'?>
	*value_real3 = *(real3*)partial_alpha_l;
]], self:getEnv()),
		},

		{
			name = 'DBar2_alpha_ll',
			code = template([[
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial2'alpha'?>



<?=eqn:makePartial'epsilon_LL'?>
	_3sym3 partial_gammaHat_lll = calc_partial_gammaHat_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>


	sym3 gammaBar_uu = calc_gammaBar_uu(U, x);
<?=eqn:getCode_connBar_ull()?>

	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	sym3 DBar2_alpha_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_alpha_ll.<?=xij?> = partial2_alpha_ll[<?=ij-1?>]
<?	for k,xk in ipairs(xNames) do
?>		- partial_alpha_l[<?=k-1?>] * connBar_ull.<?=xk?>.<?=xij?>
<?	end
?>	;
<? end
?>

	*value = sym3_dot(gammaBar_uu, DBar2_alpha_ll);
]], self:getEnv()),
		},
	
		{
			name = 'tracelessPart_LL',
			type = 'sym3',
			code = template([[
	
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial2'alpha'?>

<?=eqn:makePartial'epsilon_LL'?>
	_3sym3 partial_gammaHat_lll = calc_partial_gammaHat_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>

	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);

<?=eqn:getCode_connBar_ull()?>

	//DBar_i DBar_j alpha = alpha,ij - connBar^k_ij alpha,k
	sym3 DBar2_alpha_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_alpha_ll.<?=xij?> = partial2_alpha_ll[<?=ij-1?>]
<?	for k,xk in ipairs(xNames) do
?>		- partial_alpha_l[<?=k-1?>] * connBar_ull.<?=xk?>.<?=xij?>
<?	end
?>	;
<? end
?>

<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
	
	real3 partial_phi_l;
	sym3 partial2_phi_ll;
	{
		
		//partial_phi_l.i := phi_,i = -W_,i / (2 W) 
<? for i,xi in ipairs(xNames) do
?>		partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
<? end ?>

		//This is only used by ABar_ij,t:
		//partial2_phi_ll.ij := phi_,ij = 1/(2W) (-W_,ij + W_,i W_,j / W)
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
?>		partial2_phi_ll.<?=xij?> = .5 * (
				-partial2_W_ll[<?=ij-1?>] 
				+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
			) / U->W;
<? end ?>
	}

	
	//DBar2_phi_ll.ij := DBar_i DBar_j phi = phi_,ij - connBar^k_ij phi_,k
	sym3 DBar2_phi_ll;
<? for ij,xij in ipairs(symNames) do
?>	DBar2_phi_ll.<?=xij?> = partial2_phi_ll.<?=xij?> 
<?	for k,xk in ipairs(xNames) do	
?>		- connBar_ull.<?=xk?>.<?=xij?> * partial_phi_l.<?=xk?>
<?	end
?>	;
<? end
?>

	
	sym3 tracelessPart_ll;
<? for ij,xij in ipairs(symNames) do
	local i,j = from6to3x3(ij)
	local xi,xj = xNames[i],xNames[j]
?>	tracelessPart_ll.<?=xij?> = 0.
			+ 2. * partial_phi_l.<?=xi?> * partial_alpha_l[<?=j-1?>]
			+ 2. * partial_phi_l.<?=xj?> * partial_alpha_l[<?=i-1?>]
			+ U->alpha * (0.
				- 2. * DBar2_phi_ll.<?=xij?>
				+ 4. * partial_phi_l.<?=xi?> * partial_phi_l.<?=xj?>
				- 8. * M_PI * U->S_ll.<?=xij?>
			)
		;
<? end
?>
	tracelessPart_ll = tracefree(tracelessPart_ll, gammaBar_ll, gammaBar_uu);

	*value_sym3 = tracelessPart_ll; 
]], self:getEnv()),
		},
--]=]	
	}

	local derivOrder = 2 * self.solver.numGhost
--[=[
--[[ expansion:
2003 Thornburg:  ... from Wald ...
Theta = n^i_;i + K_ij n^i n^j - K
= n^i_,i + Gamma^i_ji n^j + K_ij (n^i n^j - gamma^ij)
... in ADM: n^i = -beta^i / alpha ...
= (-beta^i / alpha)_,i + Gamma^i_ji (-beta^j / alpha) + K_ij (beta^i beta^j / alpha^2 - gamma^ij)
= -beta^i_,i / alpha
	+ beta^i alpha_,i / alpha^2
	- beta^i (1/2 |g|_,i / |g|) / alpha
	+ K_ij beta^i beta^j / alpha^2
	- K

Gamma^j_ij = (ln sqrt(gamma))_,i 
= .5 (ln gamma)_,i 
= .5 gamma_,i / gamma
using gamma = gammaHat / W^6
= .5 (gammaHat W^-6)_,i / (gammaHat W^-6)
= .5 (gammaHat_,i W^-6 - 3 gammaHat W^-7) / (gammaHat W^-6)
= .5 gammaHat_,i / gammaHat - 3 / W
= GammaHat^j_ij - 3 / W
--]]
	vars:insert{
		name = 'expansion', 
		code = template([[
<?=eqn:makePartial'W'?>
<?=eqn:makePartial'alpha'?>
<?=eqn:makePartial'beta_U'?>
for i,xi in ipairs(xNames) do
?> + partial_beta_ul[<?=i-1?>].<?=xi?><?
end ?>;

	real exp_4phi = 1. / calc_exp_neg4phi(U);

	//gamma_ij = exp(4 phi) gammaBar_ij
	sym3 gamma_ll = sym3_real_mul(calc_gammaBar_ll(U, x), exp_4phi);

	//K_ij = exp(4 phi) ABar_ij + 1/3 gamma_ij K 
	sym3 K_ll = sym3_add(
		sym3_real_mul(U->ABar_ll, exp_4phi),
		sym3_real_mul(gamma_ll, U->K/3.));

	*value = -tr_partial_beta / U->alpha
<? 
for i,xi in ipairs(xNames) do
?>		+ U->beta_u.<?=xi?> * partial_alpha_l[<?=i-1?>] / (U->alpha * U->alpha) 
		- U->beta_u.<?=xi?> * partial_alpha_l[<?=i-1?>] / (U->alpha * U->alpha) 
		+ 3. * partial_W_l[<?=i-1?>] * U->beta_u.<?=xi?> / (U->W * U->alpha)
<?	for j,xj in ipairs(xNames) do
?>		+ K_ll.<?=sym(i,j)?> * U->beta_u.<?=xi?> * U->beta_u.<?=xj?> / (U->alpha * U->alpha)
<?	end
end
?>		- U->K;
]], self:getEnv())
	}
--]=]		
--[=[
		--[[ ADM geodesic equation spatial terms:
		-Gamma^i_tt = 
			- gamma^ij alpha_,j

			+ alpha^-1 (
				gamma^ij beta^l gamma_kl beta^k_,j
				+ 1/2 gamma^ij gamma_kl,j beta^k beta^l
				- beta^i_,t
				- gamma^ij beta^k gamma_jk,t

				+ alpha^-1 beta^i (
					alpha_,t
					+ beta^j alpha_,j

					+ alpha^-1 (
						beta^i 1/2 beta^j beta^k gamma_jk,t
						- beta^i 1/2 beta^j beta^k beta^l gamma_kl,j
						- beta^i beta^j beta^l gamma_kl beta^k_,j
					)
				)
			)

		substitute 
		alpha_,t = -alpha^2 f K + beta^j alpha_,j
		beta^k_,t = B^k
		gamma_jk,t = -2 alpha K_jk + gamma_jk,l beta^l + gamma_lj beta^l_,k + gamma_lk beta^l_,j
		--]]
	vars:insert{
		name = 'gravity',
		code = template([[
<?=eqn:makePartial'alpha'?>

	real _1_alpha = 1. / U->alpha;

	sym3 gamma_uu = calc_gamma_uu(U, x);
	real3 partial_alpha_u = sym3_real3_mul(gamma_uu, *(real3*)partial_alpha_l);		//alpha_,j gamma^ij = alpha^,i
	
<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> = -partial_alpha_u.<?=xi?>;
<? end
?>

<? if eqn.useShift ~= 'none' then ?>

<?=eqn:makePartial'beta_U'?>
<?=eqn:makePartial'epsilon_LL'?>
	
	//W = exp(-2 phi)
	real _1_W = 1. / U->W;
	
	//gamma_ij = W^-2 gammaBar_ij
	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	sym3 gamma_ll = sym3_real_mul(gammaBar_ll, _1_W * _1_W);
	
	//gamma_ij,k = W^-2 gammaBar_ij,k - 2 W^-3 gammaBar_ij W_,k
	<?=eqn:makePartial'W'?>
	_3sym3 partial_gamma_lll = {
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = sym3_sub(
			sym3_real_mul(partial_epsilon_lll[<?=i-1?>], _1_W * _1_W),
			sym3_real_mul(gammaBar_ll, 2. * partial_W_l[<?=i-1?>] * _1_W * _1_W * _1_W)),
<? end
?>	};

	//TODO
	real dt_alpha = 0.;
	sym3 dt_gamma_ll = sym3_zero;
	
	real partial_alpha_dot_beta = real3_dot(U->beta_u, *(real3*)partial_alpha_l);	//beta^j alpha_,j

	real3 beta_l = sym3_real3_mul(gamma_ll, U->beta_u);								//beta^j gamma_ij
	real3 beta_dt_gamma_l = sym3_real3_mul(dt_gamma_ll, U->beta_u);					//beta^j gamma_ij,t
	real beta_beta_dt_gamma = real3_dot(U->beta_u, beta_dt_gamma_l);				//beta^i beta^j gamma_ij,t
	
	real3 beta_dt_gamma_u = sym3_real3_mul(gamma_uu, beta_dt_gamma_l);				//gamma^ij gamma_jk,t beta^k

	//beta^i beta^j beta^k gamma_ij,k
	real beta_beta_beta_partial_gamma = 0.<?
for i,xi in ipairs(xNames) do
?> + U->beta_u.<?=xi?> * real3_weightedLenSq(U->beta_u, partial_gamma_lll.<?=xi?>)<?
end ?>;

	//beta_j beta^j_,i
	real3 beta_dbeta_l;
<? for i,xi in ipairs(xNames) do
?>	beta_dbeta_l.<?=xi?> = 0.
<?	for j,xj in ipairs(xNames) do
?>		+ beta_l.<?=xj?> * partial_beta_ul[<?=i-1?>].<?=xj?>
<?	end
?>	;
<? end
?>

	//beta_j beta^j_,i beta^i
	real beta_beta_dbeta = real3_dot(U->beta_u, beta_dbeta_l);

	//beta_j beta^j_,k gamma^ik
	real3 beta_dbeta_u = sym3_real3_mul(gamma_uu, beta_dbeta_l);

	//gamma_kl,j beta^k beta^l
	real3 beta_beta_dgamma_l = (real3){
<? for i,xi in ipairs(xNames) do
?>		.<?=xi?> = real3_weightedLenSq(U->beta_u, partial_gamma_lll.<?=xi?>),
<? end
?>	};

	real3 beta_beta_dgamma_u = sym3_real3_mul(gamma_uu, beta_beta_dgamma_l);

<? for i,xi in ipairs(xNames) do
?>	value_real3-><?=xi?> +=
		_1_alpha * (
			beta_dbeta_u.<?=xi?>
			+ .5 * beta_beta_dgamma_u.<?=xi?>	
			- U->B_U.<?=xi?>
			- beta_dt_gamma_u.<?=xi?>

			+ _1_alpha * U->beta_u.<?=xi?> * (
				.5 * dt_alpha
				+ partial_alpha_dot_beta

				+ _1_alpha * (
					.5 * beta_beta_dt_gamma
					- .5 * beta_beta_beta_partial_gamma 
					- beta_beta_dbeta
				)
			)
		)
	; 
<? end
?>
<? end	-- eqn.useShift ?>
]],			applyCommon{
				eqn = self,
				solver = self.solver,
			}
		), 
		type = 'real3',
	}
--]=]
-- [=[
	do
		vars:insert{
			name = 'RBar_LL',
			type = 'sym3',
			code = template([[
<?=eqn:makePartial'LambdaBar_U'?>
<?=eqn:makePartial'epsilon_LL'?>
<?=eqn:makePartial2'epsilon_LL'?>

	sym3 gammaBar_LL = calc_gammaBar_LL(U, x);
	real det_gammaBarLL = calc_det_gammaBarLL(x);
	sym3 gammaBar_UU = sym3_inv(gammaBar_LL, det_gammaBarLL);
	sym3 gammaBar_ll = sym3_rescaleToCoord_LL(gammaBar_LL, x);
	sym3 gammaBar_uu = sym3_rescaleToCoord_UU(gammaBar_UU, x);

	_3sym3 connHat_lll = calc_connHat_lll(x);
	_3sym3 connHat_ull = calc_connHat_ull(x);

	_3sym3 partial_gammaHat_lll = calc_partial_gammaHat_lll(x);
	
	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
	// = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end
?>

<?=eqn:getCode_connBar_ull()?>
	_3sym3 connBar_ULL = _3sym3_rescaleFromCoord_ull(connBar_ull, x);

	real3 LambdaBar_u = real3_rescaleToCoord_U(U->LambdaBar_U, x);
	real3 Delta_U = real3_sub(U->LambdaBar_U, mystery_C_U);

	sym3sym3 partial2_gammaHat_llll = calc_partial2_gammaHat_llll(x);
	sym3 gammaHat_uu = calc_gammaHat_uu(x);
	
	_3sym3 partial_connHat_ulll[3];
	calc_partial_connHat_ulll(partial_connHat_ulll, x);
	
	_3sym3 connHat_LLL = _3sym3_rescaleFromCoord_lll(connHat_lll, x);
	_3sym3 connHat_ULL = _3sym3_rescaleFromCoord_ull(connHat_ull, x);
	
	_3sym3 Delta_ULL = _3sym3_sub(connBar_ULL, connHat_ULL);
	_3sym3 Delta_LLL = sym3_3sym3_mul(gammaBar_LL, Delta_ULL);

	sym3 RBar_LL;
<?=eqn:getCode_RBar_LL()?>
	*value_sym3 = RBar_LL;
]], self:getEnv()),
		}
	end
--]=]		
--[=[		
	vars:insert{
		name = 'RPhi_LL',
		type = 'sym3',
		code = template([[

<?=eqn:makePartial'W'?>
<?=eqn:makePartial2'W'?>
-
-	real3 partial_phi_l;
-<? for i,xi in ipairs(xNames) do
-?>	partial_phi_l.<?=xi?> = -partial_W_l[<?=i-1?>] / (2. * U->W);
-<? end ?>
-
-	sym3 partial2_phi_ll;
-<? for ij,xij in ipairs(symNames) do
-	local i,j = from6to3x3(ij)
-?>	partial2_phi_ll.<?=xij?> = .5 * (
-			-partial2_W_ll[<?=ij-1?>] 
-			+ partial_W_l[<?=i-1?>] * partial_W_l[<?=j-1?>] / U->W
-		) / U->W;
-<? end ?>
-	
-	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
-	real det_gammaBar = calc_det_gammaBar(x);
-	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
-
-<?=eqn:makePartial'epsilon_LL'?>
-	_3sym3 partial_gammaHat_lll = calc_partial_gammaHat_lll(x);
 	
-	//partial_gammaBar_lll.k.ij := gammaBar_ij,k
-	// = gammaHat_ij,k + epsilon_ij,k
-	_3sym3 partial_gammaBar_lll;
-<? 
-for k,xk in ipairs(xNames) do
-	for ij,xij in ipairs(symNames) do
-?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
-<?	end
-end
-?>

-<?=eqn:getCode_connBar_ull()?>

-	sym3 DBar2_phi_ll = sym3_sub(
-		partial2_phi_ll,
-		real3_3sym3_dot1(
-			partial_phi_l,
-			connBar_ull
-		)
-	);
-
-	real tr_DBar2_phi = sym3_dot(gammaBar_uu, DBar2_phi_ll);
-
-	//2008 Alcubierre eqn 2.8.18
-	//2010 Baumgarte, Shapiro eqn 3.10
-	sym3 RPhi_ll = {
-<? for ij,xij in ipairs(symNames) do
-	local i,j = from6to3x3(ij)
-	local xi,xj = xNames[i],xNames[j]
-?>		.<?=xij?> = 
-			- 2. * DBar2_phi_ll.<?=xij?>
-			+ 4. * partial_phi_l.<?=xi?> * partial_phi_l.<?=xj?>
-			+ gammaBar_ll.<?=xij?> * (
-				- 2. * tr_DBar2_phi
-				- 4. * real3_weightedLenSq(
-					partial_phi_l,
-					gammaBar_uu
-				)
-			),
-<? end
-?>	};

	*value_sym3 = RPhi_ll;
]], self:getEnv()),
	}

--[[
gammaBar_ij = gammaHat_ij + epsilon_ij
= gammaHat_ij + epsilon_IJ e_i^I e_j^J

gammaBar_ij,kl = gammaHat_ij,kl + (epsilon_IJ e_i^I e_j^J)_,kl
= gammaHat_ij,kl + (epsilon_IJ,k e_ij^IJ + epsilon_IJ e_ij^IJ_,k)_,l
= gammaHat_ij,kl
	+ epsilon_IJ,kl e_ij^IJ 
	+ epsilon_IJ,l e_ij^IJ_,k
	+ epsilon_IJ,k e_ij^IJ_,l
	+ epsilon_IJ e_ij^IJ_,kl

gammaBar^kl gammaBar_ij,kl

gammaBar^kl = inv(gammaBar_kl)
= inv(gammaHat_kl + epsilon_kl)
--]]
--]=]	

--[=[
	do
		vars:insert{
			name='del gammaBar_ll',
			type = 'sym3',
			code = template([[

	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);

<?=eqn:makePartial2'epsilon_LL'?>
	
	sym3sym3 partial2_gammaHat_llll = calc_partial2_gammaHat_llll(x);
	
	//partial2_gammaBar_llll.kl.ij := gammaBar_ij,kl
	sym3sym3 partial2_gammaBar_llll;
<? 
for ij,xij in ipairs(symNames) do
	for kl,xkl in ipairs(symNames) do
?>	partial2_gammaBar_llll.<?=xkl?>.<?=xij?> = partial2_epsilon_llll[<?=kl-1?>].<?=xij?> + partial2_gammaHat_llll.<?=xkl?>.<?=xij?>;
<?	end
end 
?>

	sym3 trBar_partial2_gammaBar_ll;
<? for ij,xij in ipairs(symNames) do
?>	trBar_partial2_gammaBar_ll.<?=xij?> = 0.
<?	for k,xk in ipairs(xNames) do
		for l,xl in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * partial2_gammaBar_llll.<?=sym(k,l)?>.<?=xij?>
<?		end
	end
?>	;
<? end
?>
	*value_sym3 = trBar_partial2_gammaBar_ll;
]], env),
		}
	end
--]=]

--[=[
	vars:insert{
		name = 'tr34 (gamma*dGamma)',
		type = 'real3x3',
		code = template([[

	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);
		
	_3sym3 partial_connHat_ulll[3];
	calc_partial_connHat_ulll(partial_connHat_ulll, x);

	real3x3 tr34_gamma_dGamma_ll;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	tr34_gamma_dGamma_ll.<?=xi?>.<?=xj?> = 0.
<?		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
				for m,xm in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * gammaBar_ll.<?=sym(i,m)?> * partial_connHat_ulll[<?=l-1?>].<?=xm?>.<?=sym(j,k)?>
<?				end
			end
		end
?>	;
<?	end
end
?>
	
	*value_real3x3 = tr34_gamma_dGamma_ll;
]], self:getEnv()),
	}
--]=]

--[=[
	vars:insert{
		name = 'tr14 (Gamma*dgamma)',
		type = 'real3x3',
		code = template([[
<?=eqn:makePartial'epsilon_LL'?>

	_3sym3 partial_gammaHat_lll = calc_partial_gammaHat_lll(x);

	sym3 gammaBar_ll = calc_gammaBar_ll(U, x);
	real det_gammaBar = calc_det_gammaBar(x);
	sym3 gammaBar_uu = sym3_inv(gammaBar_ll, det_gammaBar);

	_3sym3 connHat_ull = calc_connHat_ull(x);

	//partial_gammaBar_lll.k.ij := gammaBar_ij,k = gammaHat_ij,k + epsilon_ij,k
	_3sym3 partial_gammaBar_lll;
<? 
for k,xk in ipairs(xNames) do
	for ij,xij in ipairs(symNames) do
?>	partial_gammaBar_lll.<?=xk?>.<?=xij?> = partial_epsilon_lll[<?=k-1?>].<?=xij?> + partial_gammaHat_lll.<?=xk?>.<?=xij?>;
<?	end
end

	real3x3 tr14_Gamma_dgamma_ll;
<? 
for i,xi in ipairs(xNames) do
	for j,xj in ipairs(xNames) do
?>	tr14_Gamma_dgamma_ll.<?=xi?>.<?=xj?> = 0.
<?		for k,xk in ipairs(xNames) do
			for l,xl in ipairs(xNames) do
				for m,xm in ipairs(xNames) do
?>		+ gammaBar_uu.<?=sym(k,l)?> * partial_gammaBar_lll.<?=xl?>.<?=sym(m,i)?> * connHat_ull.<?=xm?>.<?=sym(k,j)?>
<?				end
			end
		end
?>	;
<?	end
end
?>
	*value_real3x3 = tr14_Gamma_dgamma_ll;
]], self:getEnv()),
	}
--]=]

	--[[ hmm? not working.
	vars:insert{name='x', type='real3', code='*value_real3=x;'}
	--]]
	-- [[
	vars:insert{name='x', code='*value=x.x;'}
	vars:insert{name='y', code='*value=x.y;'}
	vars:insert{name='z', code='*value=x.z;'}
	--]]

	return vars
end

function BSSNOKFiniteDifferenceEquation:fillRandom(epsilon)
	local ptr = BSSNOKFiniteDifferenceEquation.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.numCells-1 do
		ptr[i].alpha = ptr[i].alpha + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return BSSNOKFiniteDifferenceEquation

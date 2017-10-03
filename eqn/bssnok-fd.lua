local file = require 'ext.file'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local symmath = require 'symmath'
local clnumber = require 'cl.obj.number'
local NumRelEqn = require 'eqn.numrel'
local makestruct = require 'eqn.makestruct'
require 'common'(_G)

local typeInfo = {
	real = {
		add = function(a,b) 
			if tonumber(a) == 0 then return b end
			if tonumber(b) == 0 then return a end
			return '('..a..') + ('..b..')' 
		end, 
		sub = function(a,b) return '('..a..') - ('..b..')' end, 
		scale = function(a,b) 
			if tonumber(a) == 1 then return b end
			if tonumber(b) == 1 then return a end
			return '('..a..') * ('..b..')' 
		end, 
		zero = '0',
	},
	real3 = {
		add = function(a,b) return 'real3_add('..a..', '..b..')' end,
		sub = function(a,b) return 'real3_sub('..a..', '..b..')' end,
		scale = function(a,b) return 'real3_scale('..a..', '..b..')' end,
		zero = '_real3(0,0,0)',
	},
	sym3 = {
		add = function(a,b) return 'sym3_add('..a..', '..b..')' end,
		sub = function(a,b) return 'sym3_sub('..a..', '..b..')' end,
		scale = function(a,b) return 'sym3_scale('..a..', '..b..')' end,
		zero = '_sym3(0,0,0,0,0,0)',
	},
}

-- derivCoeffs[derivative][accuracy] = {coeffs...}
local derivCoeffs = {
	-- antisymmetric coefficients 
	{
		[2] = {.5},
		[4] = {2/3, -1/12},
		[6] = {3/4, -3/20, 1/60},
		[8] = {4/5, -1/5, 4/105, -1/280},
	},
	-- symmetric
	{
		[2] = {[0] = -2, 1},
		[4] = {[0] = -5/2, 4/3, -1/12},
		[6] = {[0] = -49/18, 3/2, -3/20, 1/90},
		[8] = {[0] = -205/72, 8/5, -1/5, 8/315, -1/560},
	},
}

local function makePartial(order, solver, field, fieldType)
	local suffix = 'l'
	if not field:find'_' then suffix = '_' .. suffix end
	local name = 'partial_'..field..suffix
	local fieldTypeInfo = assert(typeInfo[fieldType], "failed to find typeInfo for "..fieldType)
	local add, sub, scale, zero = fieldTypeInfo.add, fieldTypeInfo.sub, fieldTypeInfo.scale, fieldTypeInfo.zero
	local d1coeffs = assert(derivCoeffs[1][order], "couldn't find 1st derivative coefficients of order "..order)
	local lines = table{'\t'..fieldType..' '..name..'[3];\n'}
	for i,xi in ipairs(xNames) do
		local namei = name..'['..(i-1)..']'
		local expr = zero
		if i <= solver.dim then
			for j,coeff in ipairs(d1coeffs) do
				expr = add(expr, scale(sub(
						'U['..j..' * stepsize.'..xi..'].'..field,
						'U[-'..j..' * stepsize.'..xi..'].'..field
					), clnumber(coeff)))
			end
			expr = scale(expr, '1. / grid_dx'..(i-1))
		end
		lines:insert('\t'..namei..' = '..expr..';')
	end
	return lines:concat'\n'
end

local function makePartial2(order, solver, field, fieldType, nameOverride)
	local suffix = 'll'
	if not field:find'_' then suffix = '_' .. suffix end
	local name = nameOverride or ('partial2_'..field..suffix)
	local fieldTypeInfo = assert(typeInfo[fieldType], "failed to find typeInfo for "..fieldType)
	local add, sub, scale, zero = fieldTypeInfo.add, fieldTypeInfo.sub, fieldTypeInfo.scale, fieldTypeInfo.zero
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
			local expr = scale('U->'..field, d2coeffs[0])
			for k,coeff in ipairs(d2coeffs) do
				expr = add(
					expr, 
					scale(
						add(
							'U['..k..' * stepsize.s'..(i-1)..'].'..field,
							'U[-'..k..' * stepsize.s'..(i-1)..'].'..field),
						clnumber(coeff)))
			end
			expr = scale(expr, '1. / (grid_dx'..(i-1)..' * grid_dx'..(i-1)..')')
			lines:insert('\t'..nameij..' = '..expr..';')
		else
			local expr = zero
			for k,coeff_k in ipairs(d1coeffs) do
				for l,coeff_l in ipairs(d1coeffs) do
					expr = add(expr, scale(
						sub(
							add(
								'U['..k..' * stepsize.'..xi..' + '..l..' * stepsize.'..xj..'].'..field,
								'U[-'..k..' * stepsize.'..xi..' - '..l..' * stepsize.'..xj..'].'..field),
							add(
								'U[-'..k..' * stepsize.'..xi..' + '..l..' * stepsize.'..xi..'].'..field,
								'U['..k..' * stepsize.'..xi..' - '..l..' * stepsize.'..xi..'].'..field)), 
						clnumber(coeff_k * coeff_l)))
				end
			end
			expr = scale(expr, '1. / (grid_dx'..(i-1)..' * grid_dx'..(i-1)..')')
			lines:insert('\t'..nameij..' = '..expr..';')
		end
	end
	return lines:concat'\n'
end


local BSSNOKFiniteDifferenceEquation = class(NumRelEqn)

BSSNOKFiniteDifferenceEquation.name = 'BSSNOK finite difference' 

BSSNOKFiniteDifferenceEquation.hasEigenCode = true

-- options:

-- needs to be defined up front
-- otherwise rebuild intVars based on it ...
BSSNOKFiniteDifferenceEquation.useHypGammaDriver = true

-- use chi = 1/psi instead of phi, as described in 2006 Campanelli 
-- it should be used with the hyperbolic gamma driver
BSSNOKFiniteDifferenceEquation.useChi = true

local intVars = table{
	{alpha = 'real'},			-- 1
	{beta_u = 'real3'},         -- 3: beta^i
	{gammaBar_ll = 'sym3'},    -- 6: gammaBar_ij, only 5 dof since det gammaBar_ij = 1
                                                                                                 
	BSSNOKFiniteDifferenceEquation.useChi 
		and {chi = 'real'}		-- 1
		or {phi = 'real'},		-- 1
	
	{K = 'real'},               -- 1
	{ATilde_ll = 'sym3'},       -- 6: ATilde_ij, only 5 dof since ATilde^k_k = 0
	{connBar_u = 'real3'},      -- 3: connBar^i = gammaBar^jk connBar^i_jk = -partial_j gammaBar^ij
}

if BSSNOKFiniteDifferenceEquation.useHypGammaDriver then
	intVars:insert{B_u = 'real3'}
end

local consVars = table()
:append(intVars)
:append{
	--hyperbolic variables:
	--real3 a;			//3: a_i
	--sym3 dTilde[3];		//18: dTilde_ijk, only 15 dof since dTilde_ij^j = 0
	--real3 Phi;			//3: Phi_i

	--stress-energy variables:
	{rho = 'real'},		--1: n_a n_b T^ab
	{S_u = 'real3'},			--3: -gamma^ij n_a T_aj
	{S_ll = 'sym3'},			--6: gamma_i^c gamma_j^d T_cd

	--constraints:
	{H = 'real'},				--1
	{M_u = 'real3'},			--3

	-- aux variable
	{gammaBar_uu = 'sym3'},		--6
}

BSSNOKFiniteDifferenceEquation.consVars = consVars
BSSNOKFiniteDifferenceEquation.numIntStates = makestruct.countReals(intVars)

BSSNOKFiniteDifferenceEquation.useConstrainU = true
BSSNOKFiniteDifferenceEquation.useSourceTerm = true


function BSSNOKFiniteDifferenceEquation:createInitState()
	BSSNOKFiniteDifferenceEquation.super.createInitState(self)
	self:addGuiVars{
		{name='constrain_det_gammaBar_ll', value=true},
		{name='constrain_tr_ATilde_ll', value=true},
		{name='useGammaDriver', value=false},
		{name='diffuseSigma', value=.01},
	}
end

function BSSNOKFiniteDifferenceEquation:getTemplateEnv()
	local derivOrder = 2 * self.solver.numGhost
	return {
		clnumber = clnumber,
		eqn = self,
		solver = self.solver,
		xNames = xNames,
		symNames = symNames,
		from3x3to6 = from3x3to6,
		from6to3x3 = from6to3x3,
		sym = sym,
		typeInfo = typeInfo,
		makePartial = function(...) return makePartial(derivOrder, self.solver, ...) end,
		makePartial2 = function(...) return makePartial2(derivOrder, self.solver, ...) end,
	}
end

-- should this be getInitStateCode like in eqn/euler?
function BSSNOKFiniteDifferenceEquation:getCodePrefix()
	local lines = table()
	
	lines:insert(BSSNOKFiniteDifferenceEquation.super.getCodePrefix(self))
	
	lines:insert(template([[
void setFlatSpace(global <?=eqn.cons_t?>* U) {
	U->alpha = 1.;
	U->beta_u = _real3(0,0,0);
	U->gammaBar_ll = _sym3(1,0,0,1,0,1);
<? if eqn.useChi then 
?>	U->chi = 1;
<? else
?>	U->phi = 0;
<? end
?>	U->K = 0;
	U->ATilde_ll = _sym3(1,0,0,1,0,1);
	U->connBar_u = _real3(0,0,0);
<? if eqn.useHypGammaDriver then
?>	U->B_u = _real3(0,0,0);
<? end
?>	U->gammaBar_uu = _sym3(1,0,0,1,0,1);

	//what to do with the constraint vars and the source vars?
	U->rho = 0;
	U->S_u = _real3(0,0,0);
	U->S_ll = _sym3(0,0,0,0,0,0);
	U->H = 0;
	U->M_u = _real3(0,0,0);
}

<? if eqn.useChi then
?>#define calc_exp_neg4phi(U) ((U)->chi)
<? else
?>#define calc_exp_neg4phi(U) (exp(-4. * (U)->phi))
<? end
?>

real calc_det_gamma(global const <?=eqn.cons_t?>* U) {
	real exp_neg4phi = calc_exp_neg4phi(U);
	real det_gamma = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);
	return det_gamma;
}

sym3 calc_gamma_uu(global const <?=eqn.cons_t?>* U) {
	real exp_neg4phi = calc_exp_neg4phi(U);
	sym3 gamma_uu = sym3_scale(U->gammaBar_uu, exp_neg4phi);
	return gamma_uu;
}

]], {eqn=self}))
	
	return lines:concat()
end

function BSSNOKFiniteDifferenceEquation:getInitStateCode()
	return template([[
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	real3 mids = real3_scale(real3_add(mins, maxs), .5);
	
	global <?=eqn.cons_t?>* U = UBuf + index;

	real alpha = 1.;
	real3 beta_u = _real3(0,0,0);
	sym3 gamma_ll = _sym3(1,0,0,1,0,1);
	sym3 K_ll = _sym3(0,0,0,0,0,0);
	real rho = 0.;

	<?=code?>

	U->alpha = alpha;
	U->beta_u = beta_u;

	real det_gamma = sym3_det(gamma_ll);
	sym3 gamma_uu = sym3_inv(gamma_ll, det_gamma);

	//gammaBar_ij = e^(-4phi) gamma_ij
	//real exp_neg4phi = exp(-4 * U->phi);
	real exp_neg4phi = 1./cbrt(det_gamma);

<? if eqn.useChi then 
?>	U->chi = exp_neg4phi;
<? else
?>	U->phi = log(det_gamma) / 12.;
<? end 
?>
	U->gammaBar_ll = sym3_scale(gamma_ll, exp_neg4phi);
	U->gammaBar_uu = sym3_inv(U->gammaBar_ll, 1.);

]]--[[
<? for _,x in ipairs(xNames) do
?>	U->a.<?=x?> = calc_a_<?=x?>(x.x, x.y, x.z);
<? end
?>	
]]..[[	

	U->K = sym3_dot(K_ll, gamma_uu);
	sym3 A_ll = sym3_sub(K_ll, sym3_scale(gamma_ll, 1./3. * U->K));
	U->ATilde_ll = sym3_scale(A_ll, exp_neg4phi);
	
	U->rho = rho;
	U->S_u = _real3(0,0,0);
	U->S_ll = _sym3(0,0,0,0,0,0);
	
	U->H = 0.;
	U->M_u = _real3(0,0,0);
}

//after popularing gammaBar_ll, use its finite-difference derivative to initialize connBar_u
kernel void initDerivs(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(numGhost,numGhost);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	const global <?=eqn.cons_t?>* Up[dim];
	const global <?=eqn.cons_t?>* Um[dim];
	<? for j=0,solver.dim-1 do ?>{
		Up[<?=j?>] = U + stepsize.s<?=j?>;
		Um[<?=j?>] = U - stepsize.s<?=j?>;
	}<? end ?>

<?=makePartial('gammaBar_uu', 'sym3')?>

	//connBar^i = -gammaBar^ij_,j
<? for i,xi in ipairs(xNames) do
?>	U->connBar_u.<?=xi?> =<?
	for j,xj in ipairs(xNames) do
?> - partial_gammaBar_uul[<?=j-1?>].<?=sym(i,j)?><?
	end
?>;
<? end
?>
}
]], table(self:getTemplateEnv(), {
	code = self.initState.initState 
		and self.initState:initState(self.solver) 
		or '//no code from InitCond:initState() was provided',
}))
end

function BSSNOKFiniteDifferenceEquation:getSolverCode()
	return template(file['eqn/bssnok-fd.cl'], self:getTemplateEnv())
end

function BSSNOKFiniteDifferenceEquation:getDisplayVarCodePrefix()
	return template([[
	const global <?=eqn.cons_t?>* U = buf + index;
]], {
		eqn = self,
	})
end

function BSSNOKFiniteDifferenceEquation:getEigenTypeCode()
	return template([[
typedef struct { char unused; } <?=eqn.eigen_t?>;
]], {eqn=self})
end

function BSSNOKFiniteDifferenceEquation:getDisplayVars()	
	local vars = BSSNOKFiniteDifferenceEquation.super.getDisplayVars(self)

	local function addreal3(name)
		-- TODO override solver's real3 display var construction
		vars:insert{['|'..name..'| weighted'] = '*value = real3_weightedLen(U->'..name..', U->gammaBar_ll) / calc_exp_neg4phi(U);'}	
	end

	-- hmm, how to do the weighting stuff with gammaBar_ll ... 
	-- also, how to determine which metric to raise by ... gamma vs gammaBar
	local function addsym3_ll(name)
		vars:insert{['tr '..name..' weighted'] = '*value = sym3_dot(U->gammaBar_uu, U->'..name..') / calc_det_gamma(U);'}
	end

	addreal3'beta_u'
	addsym3_ll'gammaBar_ll'
	vars:insert{['det gammaBar-1'] = [[*value = -1. + sym3_det(U->gammaBar_ll);]]}	-- for logarithmic displays
	vars:insert{['det gamma based on phi'] = [[
	real exp_neg4phi = calc_exp_neg4phi(U);
	*value = 1. / (exp_neg4phi * exp_neg4phi * exp_neg4phi);   
]]}
	addsym3_ll'ATilde_ll'
	addreal3'connBar_u'
	if self.useHypGammaDriver then
		addreal3'B_u'
	end
	--[[
	addreal3'a'
	addsym3_ll'dTilde[0]'
	addsym3_ll'dTilde[1]'
	addsym3_ll'dTilde[2]'
	addreal3'Phi'
	--]]
	addreal3'S_u'
	addsym3_ll'S_ll'

	addreal3'M_u'

	local derivOrder = 2 * self.solver.numGhost
	vars:append{
		{S = '*value = sym3_dot(U->S_ll, calc_gamma_uu(U));'},
		{volume = '*value = U->alpha * calc_det_gamma(U);'},
		{expansion = '*value = -U->alpha * U->K;'},
		{f = '*value = calc_f(U->alpha);'},
		{['df/dalpha'] = '*value = calc_dalpha_f(U->alpha);'},
		{gamma_x = '*valuevec = real3_scale(sym3_x(U->gammaBar_ll), 1./calc_exp_neg4phi(U));', type='real3'},
		{gamma_y = '*valuevec = real3_scale(sym3_y(U->gammaBar_ll), 1./calc_exp_neg4phi(U));', type='real3'},
		{gamma_z = '*valuevec = real3_scale(sym3_z(U->gammaBar_ll), 1./calc_exp_neg4phi(U));', type='real3'},
		{K_x = '*valuevec = real3_add(sym3_x(U->ATilde_ll), real3_scale(sym3_x(U->gammaBar_ll), U->K/3.));', type='real3'},
		{K_y = '*valuevec = real3_add(sym3_y(U->ATilde_ll), real3_scale(sym3_y(U->gammaBar_ll), U->K/3.));', type='real3'},
		{K_z = '*valuevec = real3_add(sym3_z(U->ATilde_ll), real3_scale(sym3_z(U->gammaBar_ll), U->K/3.));', type='real3'},

		{
			gravity = template([[
<?=makePartial('alpha', 'real')?>
	*valuevec = real3_scale(sym3_real3_mul(calc_gamma_uu(U), *(real3*)partial_alpha_l), 1. / U->alpha);
]],				{
					eqn = self,
					solver = self.solver,
					makePartial = function(...) return makePartial(derivOrder, self.solver, ...) end,
				}
			), 
			type = 'real3',
		},
	}
	
	return vars
end

function BSSNOKFiniteDifferenceEquation:fillRandom(epsilon)
	local ptr = BSSNOKFiniteDifferenceEquation.super.fillRandom(self, epsilon)
	local solver = self.solver
	for i=0,solver.volume-1 do
		ptr[i].alpha = ptr[i].alpha + 1
		ptr[i].gammaBar_ll.xx = ptr[i].gammaBar_ll.xx + 1
		ptr[i].gammaBar_ll.yy = ptr[i].gammaBar_ll.yy + 1
		ptr[i].gammaBar_ll.zz = ptr[i].gammaBar_ll.zz + 1
	end
	solver.UBufObj:fromCPU(ptr)
	return ptr
end

return BSSNOKFiniteDifferenceEquation

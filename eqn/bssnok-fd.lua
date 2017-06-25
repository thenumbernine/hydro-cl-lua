local file = require 'ext.file'
local class = require 'ext.class'
local table = require 'ext.table'
local template = require 'template'
local Equation = require 'eqn.eqn'
local symmath = require 'symmath'
	
local xNames = table{'x', 'y', 'z'}
local symNames = table{'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}

local BSSNOKFiniteDifferenceEquation = class(Equation)
BSSNOKFiniteDifferenceEquation.name = 'BSSNOK finite difference' 
BSSNOKFiniteDifferenceEquation.numStates = 35

function BSSNOKFiniteDifferenceEquation:getTypeCode()
	return template([[
typedef union {
	real ptr[1];
	struct {
		real alpha;			//1
		real3 beta_u;		//3: beta^i
		sym3 gammaBar_ll;	//6: gammaBar_ij, only 5 dof since det gammaBar_ij = 1
	
		real phi;			//1
		real K;				//1
		sym3 ATilde_ll;		//6: ATilde_ij, only 5 dof since ATilde^k_k = 0
		real3 connBar_u;	//3: connBar&i = gammaBar^jk connBar^i_jk = -partial_j gammaBar^ij

]]--[[		
		//hyperbolic variables:
		real3 a;			//3: a_i
		sym3 dTilde[3];		//18: dTilde_ijk, only 15 dof since dTilde_ij^j = 0
		real3 Phi;			//3: Phi_i
]]..[[
		
		//stress-energy variables:
		real rho;			//1: n_a n_b T^ab
		real3 S_u;			//3: -gamma^ij n_a T_aj
		sym3 S_ll;			//6: gamma_i^c gamma_j^d T_cd
	
		//constraints:
		real H;
		real3 M_u;
	};
} <?=eqn.prim_t?>;
typedef <?=eqn.prim_t?> <?=eqn.cons_t?>;

]], {
		eqn = self,
	})
end

BSSNOKFiniteDifferenceEquation.guiVars = {
	require 'guivar.combo'{
		name = 'f',
		options = {'1', '.49', '.5', '1.5', '1.69', '1 + 1/alpha^2'},
		-- value?
	}
}

BSSNOKFiniteDifferenceEquation.initStates = require 'init.adm'

function BSSNOKFiniteDifferenceEquation:getCodePrefix()
	local initState = self.initStates[1+self.solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..self.solver.initStatePtr[0])	
	
	local alphaVar = symmath.var'alpha'
	
	local fGuiVar = self.guiVarsForName.f
	local fCode = fGuiVar.options[fGuiVar.value[0]+1]
	local fExpr = assert(loadstring('local alpha = ... return '..fCode))(alphaVar)
	
	self.codes = initState.init(self.solver, {
		f = fExpr,
		alphaVar = alphaVar,
	})

	return table.map(self.codes, function(code,name,t)
		return 'real calc_'..name..code, #t+1
	end):concat'\n'
end

function BSSNOKFiniteDifferenceEquation:getInitStateCode()
	return template([[
<?
local table = require 'ext.table'
local from3x3to6_table = {
	{1, 2, 3},
	{2, 4, 5},
	{3, 5, 6},
}
local function from3x3to6(i,j)
	return from3x3to6_table[i][j]
end
local from6to3x3_table = {{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}}
local function from6to3x3(i)
	return table.unpack(from6to3x3_table[i])
end

local function sym(a,b)
	assert(a >= 1 and a <= 3, "tried to index sym with "..a..", "..b)
	assert(b >= 1 and b <= 3, "tried to index sym with "..a..", "..b)
	if a > b then a,b = b,a end
	return xNames[a]..xNames[b]
end
?>

kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(2,2);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;

	U->alpha = calc_alpha(x.x, x.y, x.z);

	U->beta_u = _real3(0,0,0);

	sym3 gamma_ll = {
<? for _,xij in ipairs(symNames) do
?>		.<?=xij?> = calc_gamma_<?=xij?>(x.x, x.y, x.z),
<? end
?>	};
	real det_gamma = sym3_det(gamma_ll);
	sym3 gamma_uu = sym3_inv(det_gamma, gamma_ll);
	U->phi = log(det_gamma) / 12.;

	//gammaBar_ij = e^(-4phi) gamma_ij
	//real exp_neg4phi = exp(-4 * U->phi);
	real exp_neg4phi = 1./cbrt(det_gamma);
	U->gammaBar_ll = sym3_scale(gamma_ll, exp_neg4phi);

]]--[[
<? for _,x in ipairs(xNames) do
?>	U->a.<?=x?> = calc_a_<?=x?>(x.x, x.y, x.z);
<? end
?>	
]]..[[	

	sym3 K_ll = {
<? for _,xij in ipairs(symNames) do
?>		.<?=xij?> = calc_K_<?=xij?>(x.x, x.y, x.z),
<? end	
?>	};
	U->K = sym3_dot(K_ll, gamma_uu);
	sym3 A_ll = sym3_sub(K_ll, sym3_scale(gamma_ll, 1./3. * U->K));
	U->ATilde_ll = sym3_scale(A_ll, exp_neg4phi);
	
	U->rho = 0;
	U->S_u = _real3(0,0,0);
	U->S_ll = (sym3){.s={0,0,0,0,0,0}};
	
	U->H = 0.;
	U->M_u = _real3(0,0,0);
}

kernel void init_connBarU(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(2,2);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;
	
	const global <?=eqn.cons_t?>* Up[dim];
	const global <?=eqn.cons_t?>* Um[dim];
	for (int j = 0; j < dim; ++j) {
		Up[j] = U + stepsize[j];
		Um[j] = U - stepsize[j];
	}

	sym3 Up_gammaBar_ll, Um_gammaBar_ll;
	real det_Up_gammaBar_ll, det_Um_gammaBar_ll;
	sym3 Up_gammaBar_uu, Um_gammaBar_uu;
	sym3 partial_gammaBar_uul[3];
<? for i=1,solver.dim do
?>	
	Up_gammaBar_ll = Up[<?=i-1?>]->gammaBar_ll;
	Um_gammaBar_ll = Um[<?=i-1?>]->gammaBar_ll;
	
	//TODO hmm, det is supposed to be 1 ... but it isn't ... 
	det_Up_gammaBar_ll = sym3_det(Up_gammaBar_ll); 
	det_Um_gammaBar_ll = sym3_det(Um_gammaBar_ll); 
	
	Up_gammaBar_uu = sym3_inv(det_Up_gammaBar_ll, Up_gammaBar_ll);
	Um_gammaBar_uu = sym3_inv(det_Um_gammaBar_ll, Um_gammaBar_ll);

<? 	for jk,xjk in ipairs(symNames) do
?>	partial_gammaBar_uul[<?=i-1?>].<?=xjk?> = 
		(Up_gammaBar_uu.<?=xjk?> - Um_gammaBar_uu.<?=xjk?>)
			/ (2. * grid_dx<?=i-1?>);
<?	end 
end
for i=solver.dim+1,3 do
?>	partial_gammaBar_uul[<?=i-1?>] = (sym3){.s={0,0,0,0,0,0}};
<? end
?>

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
]], {
		eqn = self,
		solver = self.solver,
		symNames = symNames,
		xNames = xNames,
	})
end

function BSSNOKFiniteDifferenceEquation:getSolverCode()
	return template(file['eqn/bssnok-fd.cl'], {
		eqn = self,
		solver = self.solver,
		xNames = xNames,
		symNames = symNames,
	})
end

function BSSNOKFiniteDifferenceEquation:getDisplayVarCodePrefix()
	return template([[
	const global <?=eqn.cons_t?>* U = buf + index;
]], {
		eqn = self,
	})
end

function BSSNOKFiniteDifferenceEquation:getDisplayVars()
	local vars = table()

	local function addvar(name)
		vars:insert{[name] = 'value = U->'..name..';'}
	end

	local function addreal3(name)
		for _,i in ipairs(xNames) do
			addvar(name..'.'..i)
		end
	end

	local function addsym3(name)
		for _,xij in ipairs(symNames) do
			addvar(name..'.'..xij)
		end
	end

	addvar'alpha'
	addreal3'beta_u'
	addsym3'gammaBar_ll'
	addvar'phi'
	addvar'K'
	addsym3'ATilde_ll'
	addreal3'connBar_u'
	--[[
	addreal3'a'
	addsym3'dTilde[0]'
	addsym3'dTilde[1]'
	addsym3'dTilde[2]'
	addreal3'Phi'
	--]]
	addvar'rho'
	addreal3'S_u'
	addsym3'S_ll'

	addvar'H'
	addreal3'M_u'

	vars:insert{['det_gammaBar_ll_minus_1'] = [[value = -1. + sym3_det(U->gammaBar_ll);]]}
	
	vars:insert{tr_ATilde = [[
	sym3 gammaBar_uu = sym3_inv(1., U->gammaBar_ll);
	value = sym3_dot(gammaBar_uu, U->ATilde_ll);
]]}

	vars:insert{S = [[
	sym3 gammaBar_uu = sym3_inv(1., U->gammaBar_ll);
	real exp_neg4phi = exp(-4. * U->phi);
	sym3 gamma_uu = sym3_scale(gammaBar_uu, exp_neg4phi);
	value = sym3_dot(U->S_ll, gamma_uu);
]]}

	return vars
end

return BSSNOKFiniteDifferenceEquation 

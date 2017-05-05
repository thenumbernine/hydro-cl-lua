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
BSSNOKFiniteDifferenceEquation.numStates = 45

function BSSNOKFiniteDifferenceEquation:getTypeCode()
	return template([[
typedef union {
	real ptr[45];
	struct {
		real alpha;			//1
		real3 beta;			//3: beta^i
		sym3 gammaTilde;	//6: gammaTilde_ij, only 5 dof since det gammaTilde_ij = 1
	
		real phi;			//1
		real tr_K;			//1
		real3 Gamma;		//3: Gamma^i
		sym3 ATilde;		//6: ATilde_ij, only 5 dof since ATilde^k_k = 0
]]--[[		
		real3 a;			//3: a_i
		sym3 dTilde[3];		//18: dTilde_ijk, only 15 dof since dTilde_ij^j = 0
		real3 Phi;			//3: Phi_i
]]..[[	
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
kernel void initState(
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	global <?=eqn.cons_t?>* U = UBuf + index;

	U->alpha = calc_alpha(x.x, x.y, x.z);
	
	sym3 gamma = {
<? for _,xij in ipairs(symNames) do
?>		.<?=xij?> = calc_gamma_<?=xij?>(x.x, x.y, x.z),
<? end
?>	};
	real det_gamma = sym3_det(gamma);
	sym3 gammaU = sym3_inv(det_gamma, gamma);
	U->phi = log(det_gamma) / 12.;

	real exp_neg4phi = exp(-4 * U->phi);
	U->gammaTilde = sym3_scale(gamma, exp_neg4phi);
]]--[[
<? for _,x in ipairs(xNames) do
?>	U->a.<?=x?> = calc_a_<?=x?>(x.x, x.y, x.z);
<? end
?>	
]]..[[	
	sym3 K = {
<? for _,xij in ipairs(symNames) do
?>		.<?=xij?> = calc_K_<?=xij?>(x.x, x.y, x.z),
<? end	
?>	};
	U->tr_K = sym3_dot(K, gammaU);
	sym3 A = sym3_add(K, sym3_scale(gamma, -1./3. * U->tr_K));
	U->ATilde = sym3_scale(A, exp_neg4phi);
}
]], {
		eqn = self,
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
	addsym3'gammaTilde'
	addvar'phi'
	addvar'K'
	addreal3'Gamma'
	addsym3'ATilde'
	--[[
	addreal3'a'
	addsym3'dTilde[0]'
	addsym3'dTilde[1]'
	addsym3'dTilde[2]'
	addreal3'Phi'
	--]]

	return vars
end

return BSSNOKFiniteDifferenceEquation 

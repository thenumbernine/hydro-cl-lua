--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 
The first Bona-Masso formalism.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local file = require 'ext.file'
local Equation = require 'eqn.eqn'

local ADM_BonaMasso_3D = class(Equation)
ADM_BonaMasso_3D.name = 'ADM_BonaMasso_3D'

ADM_BonaMasso_3D.consVars = {
	'alpha',
	'gamma_xx', 'gamma_xy', 'gamma_xz', 'gamma_yy', 'gamma_yz', 'gamma_zz',
	'a_x', 'a_y', 'a_z',
	'd_xxx', 'd_xxy', 'd_xxz', 'd_xyy', 'd_xyz', 'd_xzz',
	'd_yxx', 'd_yxy', 'd_yxz', 'd_yyy', 'd_yyz', 'd_yzz',
	'd_zxx', 'd_zxy', 'd_zxz', 'd_zyy', 'd_zyz', 'd_zzz',
	'K_xx', 'K_xy', 'K_xz', 'K_yy', 'K_yz', 'K_zz',
	'V_x', 'V_y', 'V_z',
}
ADM_BonaMasso_3D.numStates = 7 + 30	-- should equal # consVars
ADM_BonaMasso_3D.numWaves = 30	-- skip alpha and gamma_ij

ADM_BonaMasso_3D.numStates = #ADM_BonaMasso_3D.consVars
ADM_BonaMasso_3D.displayVars = table()
	:append(ADM_BonaMasso_3D.consVars)
	:append{'volume'}

ADM_BonaMasso_3D.hasCalcDT = true
ADM_BonaMasso_3D.useSourceTerm = true

ADM_BonaMasso_3D.initStates = require 'init.adm'
ADM_BonaMasso_3D.initStateNames = table.map(ADM_BonaMasso_3D.initStates, function(state) return state.name end)

ADM_BonaMasso_3D.guiVars = table{
	require 'guivar.combo'{
		name = 'f',
		options = {'1', '1.69', '.49', '1 + 1/alpha^2'},
	}
}
ADM_BonaMasso_3D.guiVarsForName = ADM_BonaMasso_3D.guiVars:map(function(var) return var, var.name end)

function ADM_BonaMasso_3D:getTypeCode()
	return [[
typedef union {
	real s[6];
	struct {
		real xx, xy, xz, yy, yz, zz;
	};
} symmat3;

real symmat3_det(symmat3 m) {
	return m.xx * m.yy * m.zz
		+ m.xy * m.yz * m.xz
		+ m.xz * m.xy * m.yz
		- m.xz * m.yy * m.xz
		- m.yz * m.yz * m.xx
		- m.zz * m.xy * m.xy;
}

symmat3 symmat3_inv(real d, symmat3 m) {
	return (symmat3){
		.xx = (m.yy * m.zz - m.yz * m.yz) / d,
		.xy = (m.xz * m.yz - m.xy * m.zz) / d,
		.xz = (m.xy * m.yz - m.xz * m.yy) / d,
		.yy = (m.xx * m.zz - m.xz * m.xz) / d,
		.yz = (m.xz * m.xy - m.xx * m.yz) / d,
		.zz = (m.xx * m.yy - m.xy * m.xy) / d,
	};
}

typedef struct {
	real alpha;
	symmat3 gamma;
	real3 a;
	symmat3 d[3];
	symmat3 K;
	real3 V;
} cons_t;

]]
end

function ADM_BonaMasso_3D:getCodePrefix(solver)
	local initState = self.initStates[solver.initStatePtr[0]+1]
	
	local alphaVar = require 'symmath'.var'alpha'

	local fGuiVar = self.guiVarsForName.f
	local fCode = fGuiVar.options[fGuiVar.value[0]+1]
	local fExpr = assert(loadstring('local alpha = ... return '..fCode))(alphaVar)

	self.codes = initState.init(solver, {
		f = fExpr,
		alphaVar = alphaVar,
	})
	
	return table.map(self.codes, function(code,name,t)
		return 'real calc_'..name..code, #t+1
	end):concat'\n'
end
	
local xNames = table{'x', 'y', 'z'}
local symNames = table{'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}

function ADM_BonaMasso_3D:getInitStateCode(solver)
	local lines = table{
		[[
__kernel void initState(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	__global cons_t* U = UBuf + index;
]]
	}

	local function build(var)
		return '\tU->'..var..' = calc_'..var..'(x.x, x.y, x.z);'
	end

	build'alpha'
	symNames:map(function(xij) build('gamma_'..xij) end)
	xNames:map(function(xi) build('a_'..xi) end)	
	xNames:map(function(xk)
		symNames:map(function(xij) build('d_'..xk..xij) end)
	end)
	symNames:map(function(xij) build('K_'..xij) end)
	lines:insert'}'
	
	return lines:concat'\n'
end

function ADM_BonaMasso_3D:getSolverCode(solver)
	return require 'processcl'(file['eqn/adm3d.cl'], {solver=solver})
end

ADM_BonaMasso_3D.eigenVars = {'alpha', 'gammaUxx', 'gammaUxy', 'gammaUxz', 'gammaUyy', 'gammaUyz', 'gammaUzz', 'f'}
function ADM_BonaMasso_3D:getEigenInfo()
	local makeStruct = require 'eqn.makestruct'
	return {
		typeCode = [[
typedef struct {
	symmat3 gammaU;
	real f;
} eigen_t;

// I've thought of merging these two structures ... this is more proof
typedef eigen_t fluxXform_t;

]],
		code = nil,
		displayVars = {} -- working on this one
	}
end

function ADM_BonaMasso_3D:getCalcDisplayVarCode()
	return table{[[
	switch (displayVar) {
	case display_U_volume: value = U->alpha * sqrt(symmat3_det(U->gamma)); break;
]]
	}:append(table.map(self.consVars, function(var)
		local code = var:gsub('_', '.')
		local dx, rest = code:match('^d%.([xyz])([xyz][xyz])$')
		if dx then
			local i = xNames:find(dx)
			code = 'd['..(i-1)..'].'..rest
		end
		return '	case display_U_'..var..': value = U->'..code..'; break;'
	end)):append{
[[
		}
]]
	}:concat'\n'
end

return ADM_BonaMasso_3D

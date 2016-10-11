--[[
Based on Alcubierre 2008 "Introduction to 3+1 Numerical Relativity" on the chapter on hyperbolic formalisms. 
The first Bona-Masso formalism.
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

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
ADM_BonaMasso_3D.numStates = #ADM_BonaMasso_3D.consVars
ADM_BonaMasso_3D.displayVars = table()
	:append(ADM_BonaMasso_3D.consVars)
	:append{'volume'}

ADM_BonaMasso_3D.useSourceTerm = true

ADM_BonaMasso_3D.initStates = require 'init_adm'
ADM_BonaMasso_3D.initStateNames = table.map(ADM_BonaMasso_3D.initStates, function(state) return state.name end)

function ADM_BonaMasso_3D:codePrefix()
	return table.map(self.codes, function(code,name,t)
		return 'real calc_'..name..code, #t+1
	end):concat'\n'
end

ADM_BonaMasso_3D.guiVars = {'f'}
ADM_BonaMasso_3D.f = {
	value = 0,	-- 0-based index into options
	name = 'f',
	options = {'1', '1.69', '.49', '1 + 1/alpha^2'},
}

function ADM_BonaMasso_3D:getInitStateCode(solver)
	local initState = self.initStates[solver.initStatePtr[0]+1]
	
	local alphaVar = require 'symmath'.var'alpha'
	self.codes = initState.init(solver, ({
		{f = 1},
		{f = 1.69},
		{f = 1.49},
		{f = 1 + 1/alphaVar^2, alphaVar=alphaVar},
	})[self.f.value+1])

	local lines = table{
		self:codePrefix(),
		[[
__kernel void initState(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real4 x = CELL_X(i);
	__global cons_t* U = UBuf + index;
]]
	}

	local function build(var)
		return '\tU->'..var..' = calc_'..var..'(x.x, x.y, x.z);'
	end

	local xNames = table{'x', 'y', 'z'}
	local symNames = table{'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}
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

function ADM_BonaMasso_3D:solverCode()
	return table{
		self:codePrefix(),
		'#include "adm3d.cl"',
	}:concat'\n'
end



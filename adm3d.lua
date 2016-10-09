local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local BonaMassoADM3D = class(Equation)
BonaMassoADM3D.name = 'Bona-Masso ADM 3D'

BonaMassoADM3D.numStates = 37
BonaMassoADM3D.consVars = {
'alpha',
'gamma_xx', 'gamma_xy', 'gamma_xz', 'gamma_yy', 'gamma_yz', 'gamma_zz',
'a_x', 'a_y', 'a_z',
'd_xxx', 'd_xxy', 'd_xxz', 'd_xyy', 'd_xyz', 'd_xzz',
'd_yxx', 'd_yxy', 'd_yxz', 'd_yyy', 'd_yyz', 'd_yzz',
'd_zxx', 'd_zxy', 'd_zxz', 'd_zyy', 'd_zyz', 'd_zzz',
'K_xx', 'K_xy', 'K_xz', 'K_yy', 'K_yz', 'K_zz',
'V_x', 'V_y', 'V_z',
}

BonaMassoADM3D.initStates = {

}

function BonaMassoADM3D:getInitStateCode(solver)
	local symmath = require 'symmath'
	symmath.tostring = require 'symmath.tostring.SingleLine'		
	
	local x,y,z = symmath.vars('x', 'y', 'z')
	local vars = table{x,y,z}
	local alphaVar = symmath.var'alpha'

	local xc, yc, zc = 150, 150, 150
	local H = 5
	local sigma = 10
	local h = H * symmath.exp(-((x-xc)^2 + (y-yc)^2 + (z-zc)^2) / sigma^2)

	--local f = 1
	--local f = 1.69
	--local f = .49
	local f = 1 + kappa / alphaVar^2


end

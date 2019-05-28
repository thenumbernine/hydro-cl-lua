--[[
simple wave equation

d'lambertian phi = 0

phi_,tt - c^2 phi_;i^i = 0

phi_,tt - c^2 g^ij phi_;ij = 0
phi_,tt - c^2 g^ij (phi_,i)_;j = 0
phi_,tt - c^2 g^ij (phi_,ij - conn^k_ij phi_,k) = 0
phi_,tt - c^2 g^ij phi_,ij + c^2 conn^ij_j phi_,i = 0

let a = phi,t
let b_i = phi_,i 

phi_,t = a
a_,t - c^2 g^ij b_i,j = -c^2 conn^ij_j b_i
b_i,t - a,i = 0

produces:
[ a ]      [     0       -c^2 g^jk ] [ a ]      [ -c^2 conn^ij_j b_i ]
[b_i]_,t + [ -delta_i^k      0     ] [b_j]_,k = [        0           ]

flux jacobian in x dir:
[  0, -c^2 g^xx, -c^2 g^xy, -c^2 g^xz ]
[ -1,     0    ,     0    ,     0     ]
[  0,     0    ,     0    ,     0     ]
[  0,     0    ,     0    ,     0     ]

flux jacobian in y dir:
[  0, -c^2 g^xy, -c^2 g^yy, -c^2 g^yz ]
[  0,     0    ,     0    ,     0     ]
[ -1,     0    ,     0    ,     0     ]
[  0,     0    ,     0    ,     0     ]

flux jacobian in z dir:
[  0, -c^2 g^xz, -c^2 g^yz, -c^2 g^zz ]
[  0,     0    ,     0    ,     0     ]
[  0,     0    ,     0    ,     0     ]
[ -1,     0    ,     0    ,     0     ]

x dir has eigensystem:
lambda = {0, 0, -c sqrt(g^xx), c sqrt(g^xx)} 
eigenvectors = {
	{0, -g^xy, g^xx, 0},
	{0, -g^xz, 0, g^xx},
	{c sqrt(g^xx), 0, 0, 0},
	{-c sqrt(g^xx), 0, 0, 0},
}

y dir has eigensystem:
lambda = {0, 0, -c sqrt(g^yy), c sqrt(g^yy)}
eigenvectors = {
	{0, -g^yy, g^xy, 0},
	{0, -g^yz, 0, g^xy},
	{sqrt(g^yy), 0, 1, 0},
	{-sqrt(g^yy), 0, 1, 0},
}

z dir has eigensystem:
lambda = {0, 0, -c sqrt(g^zz), c sqrt(g^zz)},
eigenvectors = {
	{0, -g^yz, g^xz, 0},
	{0, -g^zz, 0, g^xz},
	{sqrt(g^zz), 0, 0, 1},
	{-sqrt(g^zz), 0, 0, 1},
}
--]]
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local template = require 'template'
local Equation = require 'eqn.eqn'

local Wave = class(Equation)
Wave.name = 'wave'

Wave.numStates = 4

Wave.mirrorVars = {{'v.x'}, {'v.y'}, {'v.z'}}

Wave.hasEigenCode = true
Wave.hasFluxFromConsCode = true
Wave.roeUseFluxFromCons = true

Wave.useSourceTerm = false

Wave.initStates = require 'init.euler'	 -- use rho as our initial condition

Wave.consVars = table{
	{phi_t = 'real'},
	{phi_i = 'real3'},
}

function Wave:createInitState()
	Wave.super.createInitState(self)
	self:addGuiVars{
		{name='wavespeed', value=1},
	}
end

Wave.initStateCode = [[
<?
local common = require 'common'()
local xNames = common.xNames
?>
kernel void initState(
	constant <?=solver.solver_t?>* solver,
	global <?=eqn.cons_t?>* UBuf
) {
	SETBOUNDS(0,0);
	real3 x = cell_x(i);
	real3 mids = real3_real_mul(real3_add(solver->mins, solver->maxs), .5);
	bool lhs = true<?
for i=1,solver.dim do
	local xi = xNames[i]
?> && x.<?=xi?> < mids.<?=xi?><?
end
?>;

	real rho = 0;
	real3 v = real3_zero;
	real P = 0;
	real3 D = real3_zero;
	real3 B = real3_zero;
	real ePot = 0;
	
	<?=code?>

	UBuf[index] = (<?=eqn.cons_t?>){
		//hmm, this might deserve its own initial conditions ...
		// for the initialization of these variables:
		//.phi_t = rho,
		.phi_t = P,
		.phi_i = real3_zero,
	};
}
]]

Wave.solverCodeFile = 'eqn/wave.cl'

Wave.eigenVars = {{unused = 'real'}}

function Wave:eigenWaveCodePrefix(side, eig, x)
	return template([[
	real c_sqrt_gU = solver->wavespeed * coord_sqrt_gU<?=side..side?>(<?=x?>);
]], {
		side = side,
		x = x,
	})
end

function Wave:consWaveCodePrefix(side, U, x, W)
	return self:eigenWaveCodePrefix(side, nil, x)
end

function Wave:consWaveCode(side, U, x, waveIndex)
	if waveIndex == 0 then
		return '-c_sqrt_gU' 
	elseif waveIndex == 1 or waveIndex == 2 then
		return '0'
	elseif waveIndex == 3 then
		return 'c_sqrt_gU' 
	end
	error'got a bad waveIndex'
end

Wave.eigenWaveCode = Wave.consWaveCode

return Wave

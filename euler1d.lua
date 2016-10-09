--[[
this has so much in common with Euler3D ...
and I don't want to keep updating the both of them ...
and I don't really care about this as much as the 3D version ...
so maybe I should have this subclass / steal from Euler3D?
--]]

local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'
local clnumber = require 'clnumber'

local Euler1D = class(Equation)
Euler1D.name = 'Euler1D'

Euler1D.numStates = 3

Euler1D.consVars = {'rho', 'mx', 'ETotal'}
Euler1D.primVars = {'rho', 'vx', 'P'}

-- notice the current cheap system just appends the dim to mirror on the var prefix
-- but in 1D, there is only 'mx', not 'my' or 'mz'
-- soo... the system will break for 2D and 3D. 
-- soo ... fix the system
Euler1D.mirrorVars = {{'mx'}, {}, {}}

Euler1D.displayVars = {
	'rho',
	'vx',
	'mx',
	'eInt',
	'eKin', 
	'eTotal', 
	'EInt', 
	'EKin', 
	'ETotal', 
	'P',
	'S', 
	'h',
	'H', 
	'hTotal',
	'HTotal',
} 

Euler1D.initStates = require 'init_euler'
Euler1D.initStateNames = table.map(Euler1D.initStates, function(info) return info.name end)

Euler1D.guiVars = {'gamma'}
Euler1D.gamma = 7/5

function Euler1D:getCodePrefix()
	return table{
		'#define gamma '..clnumber(self.gamma),
		'#define gamma_1 (gamma-1.)',
		'#define gamma_3 (gamma-3.)',
	}:concat'\n'
end

function Euler1D:getTypeCode()
	return 
		require 'makestruct'('prim_t', self.primVars) .. '\n' ..
		Euler1D.super.getTypeCode(self) 
end

function Euler1D:getInitStateCode(solver)
	local initState = self.initStates[1+solver.initStatePtr[0]]
	assert(initState, "couldn't find initState "..solver.initStatePtr[0])	
	local code = initState.init(solver)	
	
	return table{
		self:getCodePrefix(),
		[[
cons_t consFromPrim(prim_t W) {
	return (cons_t){
		.rho = W.rho,
		.mx = W.rho * W.vx,
		.ETotal = .5 * W.rho * W.vx * W.vx + W.P / gamma_1,
	};
}

__kernel void initState(
	__global cons_t* UBuf
) {
	SETBOUNDS(0,0);
	real4 x = CELL_X(i);
	real4 mids = (real).5 * (mins + maxs);
	bool lhs = x[0] < mids[0]
#if dim > 1
		&& x[1] < mids[1]
#endif
#if dim > 2
		&& x[2] < mids[2]
#endif
	;
	real rho = 0;
	real vx = 0;
	real P = 0;

]] .. code .. [[

	UBuf[index] = consFromPrim((prim_t){.rho=rho, .vx=vx, .P=P});
}
]],
	}:concat'\n'
end

function Euler1D:solverCode(solver)	
	return table{
		self:getCodePrefix(),
		'#include "euler1d.cl"',
	}:concat'\n'
end

-- TODO boundary methods, esp how to handle mirror

return Euler1D

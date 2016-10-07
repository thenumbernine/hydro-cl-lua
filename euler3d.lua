local class = require 'ext.class'
local table = require 'ext.table'
local Equation = require 'equation'

local Euler3D = class(Equation)
Euler3D.name = 'Euler3D'

Euler3D.numStates = 5

Euler3D.consVars = table{'rho', 'mx', 'my', 'mz', 'ETotal'}
Euler3D.primVars = table{'rho', 'vx', 'vy', 'vz', 'P'}
Euler3D.displayVars = table()
	:append(Euler3D.primVars)
	:append{'eInt', 'eKin', 'eTotal'} 

Euler3D.initStates = {'Sod', 'linear'}

Euler3D.gamma = 7/5

function Euler3D:getTypeCode()
	return [[

typedef union {
	struct {
		real x, y, z;
	};
	real v[3];
} real3;

typedef struct { 
	real rho;
	union {
		struct { real vx, vy, vz; };
		real3 v;
	};
	real P;
} prim_t;

typedef struct {
	real rho;
	union {
		struct { real mx, my, mz; };
		real3 m;
	};
	real ETotal;
} cons_t;

]]
end

function Euler3D:solverCode(clnumber)
	return table{
		'#define gamma '..clnumber(self.gamma),
		'#include "euler1d.cl"',
	}:concat'\n'
end

-- TODO boundary methods, esp how to handle mirror

return Euler3D

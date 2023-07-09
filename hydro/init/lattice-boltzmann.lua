local class = require 'ext.class'
local table = require 'ext.table'
local InitCond = require 'hydro.init.init'

local LatticeBoltzmannInitCond = class(InitCond)

local initConds = table{
	{
		name = 'cylinder',
		mins = {-1, -1, -1},
		maxs = {1, 1, 1},
		guiVars = {
			{name='r', value=.25},
		},
		getInitCondCode = function(self)
			local solver = assert(self.solver)
			-- TODO custom boundary.  rhs is set to zero.  lhs is U[-2] = U[2], U[-1] = U[1], and U[0] is not modified
			solver:setBoundaryMethods'freeflow'
			return [[

#define DBL_EPSILON 2.220446049250313080847e-16
#define DBL_EPS_COMP (1. - DBL_EPSILON)
solid = (real)(real3_len(xc) < initCond->r);

int ofs = 0;
int4 c = (int4)(0,0,0,0);		
for (int ofz = 0; ofz < solver->ofsmax.z; ++ofz) {
	c.z = ofz - (solver->ofsmax.z-1)/2;
	for (int ofy = 0; ofy < solver->ofsmax.y; ++ofy) {
		c.y = ofy - (solver->ofsmax.y-1)/2;
		for (int ofx = 0; ofx < solver->ofsmax.x; ++ofx) {
			c.x = ofx - (solver->ofsmax.x-1)/2;

			real const mu = 0.;
			real const sigma = 1.;
			real const randn = mu + (U->rho > .5 ? -1. : 1.) * sigma * sqrt(-log(DBL_EPS_COMP * U->nbhd[ofs]));

			nbhd[ofs] = 1. + .01 * randn;
			if (c.x == 1 && c.y == 0 && c.z == 0) {
				nbhd[ofs] += 2 * (1 + .2 * cos(2. * M_PI * xc.x));
			}

			rho += nbhd[ofs];
		}
	}
}
]]
		end,
	},
}:mapi(function(cl)
	return class(LatticeBoltzmannInitCond, cl)
end)

function LatticeBoltzmannInitCond:getList()
	return initConds
end

return LatticeBoltzmannInitCond 

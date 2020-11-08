local dim = 2
local args = {
	app = self, 
	eqn = 'euler',
	flux = 'roe',
	dim = 2,
	integrator = 'forward Euler',	
	cfl = .3/dim,
	fluxLimiter = 'superbee',
	
	coord = 'cartesian',
	coordArgs = {vectorComponent='anholonomic'},
	mins = {-1, -1, -1},
	maxs = {1, 1, 1},

	gridSize = cmdline.gridSize,
	
	boundary = {
		xmin = 'periodic',
		xmax = 'periodic',
		ymin = 'periodic',
		ymax = 'periodic',
		zmin = 'periodic',
		zmax = 'periodic',
	},
	initCond = 'shock bubble interaction',
}

self.solvers:insert(require 'hydro.solver.fvsolver'(args))

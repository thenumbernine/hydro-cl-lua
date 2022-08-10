local dim = 2
local args = {
	app = self, 
	eqn = 'wave',
	dim = dim,
	integrator = 'forward Euler',	
	cfl = .5/dim,
	
	flux = 'roe',
	fluxLimiter = cmdline.fluxLimiter,
	
	coord = cmdline.coord,
	coordArgs = {
		vectorComponent = cmdline.vectorComponent,
	},
	mins = cmdline.mins,
	maxs = cmdline.maxs,
	gridSize = cmdline.gridSize,
	boundary = cmdline.boundary,
	
	initCond = 'constant',
	initCondArgs = {v={1e-1,1e-1}},
}

self.solvers:insert(require 'hydro.solver.fvsolver'(args))

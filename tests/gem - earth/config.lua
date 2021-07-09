local dim = 3
local args = {
	app = self, 
	dim = dim,
	
	integrator = 'forward Euler',	
	cfl = .3/dim,
	fluxLimiter = 'superbee',

	initCond = 'self-gravitation - Earth',	-- validating units along with self-gravitation.
}

if cmdline.coord == 'cartesian' then
	args = table(args, {
		coord = 'cartesian',
		mins = cmdline.mins or {-1, -1, -1},
		maxs = cmdline.maxs or {1, 1, 1},
		
		gridSize = {32,32,32},
		boundary = {
			xmin = 'freeflow',
			xmax = 'freeflow',
			ymin = 'freeflow',
			ymax = 'freeflow',
			zmin = 'freeflow',
			zmax = 'freeflow',
		},
	}):setmetatable(nil)
elseif cmdline.coord == 'cylinder' then
	args = table(args, {
		coord = 'cylinder',
		--coordArgs = {vectorComponent='anholonomic'},	
		coordArgs = {vectorComponent='cartesian'},	
		mins = cmdline.mins or {0, 0, -1},
		maxs = cmdline.maxs or {1, 2*math.pi, 1},
		gridSize = {32, 32, 32},
		boundary = {
			-- r
			xmin='cylinderRMin',	-- use this when rmin=0
			xmax='freeflow',
			
			-- θ
			ymin='periodic',
			ymax='periodic',
			
			-- z
			zmin='freeflow',
			zmax='freeflow',
		},
	}):setmetatable(nil)
elseif cmdline.coord == 'sphere' then
	args = table(args, {
		coord = 'sphere',
		--coordArgs = {vectorComponent='anholonomic'},
		coordArgs = {vectorComponent='cartesian'},
		mins = cmdline.mins or {0, 0, -math.pi},
		maxs = cmdline.maxs or {8, math.pi, math.pi},
		gridSize = ({
			{160, 1, 1}, -- 1D
			{32, 32, 1}, -- 2D
			{16, 16, 16}, -- 3D
		})[dim],
		boundary = {
			xmin='sphereRMin',
			xmax='freeflow',
			ymin='sphereTheta',
			ymax='sphereTheta',
			zmin='periodic',
			zmax='periodic',
		},
	}):setmetatable(nil)
else
	error("can't handle coord=="..require 'ext.tolua'(coord))
end
self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler-lingr'})))

local args = {
	app = self, 
	dim = 3,
	
	integrator = 'forward Euler',	
	cfl = assert(tonumber(cmdline.cfl)),
	fluxLimiter = 'superbee',

	initCond = 'self-gravitation - Earth',	-- validating units along with self-gravitation.
}

if cmdline.coord == 'cartesian' then
	args = table(args, {
		coord = 'cartesian',
		mins = cmdline.mins or {-1, -1, -1},
		maxs = cmdline.maxs or {1, 1, 1},
		
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
		mins = cmdline.mins or {0, 0, -1},
		maxs = cmdline.maxs or {1, 2*math.pi, 1},
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
		mins = cmdline.mins or {0, 0, -math.pi},
		maxs = cmdline.maxs or {8, math.pi, math.pi},
		
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

args.coordArgs = args.coordArgs or {}
args.coordArgs.vectorComponent = assert(cmdline.vectorComponent)

args.gridSize = assert(cmdline.gridSize)
assert(#args.gridSize == 3)
for i=1,3 do assert(type(args.gridSize[i]) == 'number') end

self.solvers:insert(require 'hydro.solver.fvsolver'(table(args, {flux='roe', eqn='euler-lingr'})))

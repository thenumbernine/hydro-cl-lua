local table = require 'ext.table'
local configs = table{
-- [[ mesh with cartesian components
-- I'm using this as my golden standard for curvilinear grids
	{
		name = 'mesh',
		solverClassName = 'hydro.solver.meshsolver',
		solverArgs = {
			integrator = 'forward Euler',
			cfl = .25,
			flux = 'roe',
			eqn = 'euler',
			initState = 'spiral',
			-- no / donor cell flux limiter
			-- cartesian / holonomic vector components
			-- legacy to gridsolvers.  mesh doesn't need this:
			dim = 2,
			-- mesh-specific params:
			mesh = {
				type = 'polar2d',
				size = {64, 64},
			},
		},
	},
--]]
-- [[ grid with cartesian components
-- this is nearly identical to the above, which is good.	
	{
		name = 'grid-cartesian',
		solverClassName = 'hydro.solver.fvsolver',
		solverArgs = {
			integrator = 'forward Euler',
			cfl = .25,
			flux = 'roe',
			eqn = 'euler',
			initState = 'spiral',
			coord = 'cylinder',
			coordArgs = {vectorComponent = 'cartesian'},
			dim = 2,
			mins = {.1, 0, -.5},
			maxs = {1, 2*math.pi, .5},
			gridSize = {64, 64, 1},
			boundary = {
				xmin='freeflow',
				xmax='freeflow',
				ymin='periodic',
				ymax='periodic',
				zmin='freeflow',
				zmax='freeflow',
			},
		},
	},
--]]
-- [[ grid with orthonormal grid-aligend components
-- this starts to deviate ...
	{
		name = 'grid-orthonormal',
		solverClassName = 'hydro.solver.fvsolver',
		solverArgs = {
			integrator = 'forward Euler',
			cfl = .25,
			flux = 'roe',
			eqn = 'euler',
			initState = 'spiral',
			coord = 'cylinder',
			coordArgs = {vectorComponent = 'anholonomic'},
			dim = 2,
			mins = {.1, 0, -.5},
			maxs = {1, 2*math.pi, .5},
			gridSize = {64, 64, 1},
			boundary = {
				xmin='freeflow',
				xmax='freeflow',
				ymin='periodic',
				ymax='periodic',
				zmin='freeflow',
				zmax='freeflow',
			},
		},
	},
--]]
-- [[ grid with grid coordinate components
-- this encounters numericals errors
	{
		name = 'grid-coordinate',
		solverClassName = 'hydro.solver.fvsolver',
		solverArgs = {
			integrator = 'forward Euler',
			cfl = .25,
			flux = 'roe',
			eqn = 'euler',
			initState = 'spiral',
			coord = 'cylinder',
			coordArgs = {vectorComponent = 'holonomic'},
			dim = 2,
			mins = {.1, 0, -.5},
			maxs = {1, 2*math.pi, .5},
			gridSize = {64, 64, 1},
			boundary = {
				xmin='freeflow',
				xmax='freeflow',
				ymin='periodic',
				ymax='periodic',
				zmin='freeflow',
				zmax='freeflow',
			},
		},
	},
--]]
}:mapi(function(config)
	config.trackVars = {'U v'}	-- I'm not using the cmdline system ... TODO fix / merge these?
	return config
end)

return configs

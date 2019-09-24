local ffi = require 'ffi'
require 'ffi.c.stdlib'
local unistd = require 'ffi.c.unistd'

local rundirp = unistd.getcwd(nil, 0)
local rundir = ffi.string(rundirp)
ffi.C.free(rundirp)

local class = require 'ext.class'
local table = require 'ext.table'
unistd.chdir'../..'

__useConsole__ = true

local dim, solvername, initState = ...
dim = math.floor(assert(tonumber(dim)))
assert(solvername)
assert(initState)

-- 2013 Baumgarte et al, section IV A 1 example
local cfg = {
	eqn = 'bssnok-fd-num', 
	eqnArgs = {
		--useScalarField = true,	-- needed for the scalar field init cond below
	
		cflMethod = '2008 Alcubierre',
		--cflMethod = '2013 Baumgarte et al, eqn 32',
		--cflMethod = '2017 Ruchlin et al, eqn 53',
	},
	dim = dim,
	integrator = 'Runge-Kutta 4',	-- the paper says PIRK
	--integrator = 'backward Euler',
	--integratorArgs = {verbose=true},
	cfl = .1/dim,
	
	-- [[
	coord = 'sphere',
	--coord = 'sphere-log-radial',
	mins = {0, 0, -math.pi},
	maxs = {16, math.pi, math.pi},
	gridSize = assert(({
		{128, 1, 1},
		{64, 16, 1},
		{8, 8, 8},
	})[dim]),
	boundary = {
		xmin='sphereCenter',
		
		-- runs a gauge wave until t=...
		--xmax='fixed',		-- 5.3875
		--xmax='freeflow',	-- diverges near rmin after t=60 or so
		--xmax='linear',	-- 13.1875
		xmax='quadratic',	-- 10.6125
		
		ymin='spherePolar',
		ymax='spherePolar',
		zmin='periodic',
		zmax='periodic',
	},
	--]]
	--[[
	coord = 'cartesian',
	mins = {-4,-4,-4},
	maxs = {4,4,4},
	gridSize = ({
		{250, 1, 1},
		{40, 40, 1},
		{16, 16, 16},
	})[dim],

	boundary = {
		xmin='freeflow',
		xmax='freeflow',
		ymin='freeflow',
		ymax='freeflow',
		zmin='freeflow',
		zmax='freeflow',
	},
	--]]
	
	initState = initState,
}

local duration = 100

local filename = rundir..'/dim='..dim..',solver='..solvername..',init='..initState..'.txt'
local f = io.open(filename, 'w')

local App = class(require 'app')
function App:setup(clArgs)
	cfg.app = self
	local solver = require('solver.'..solvername)(cfg)
	solver.checkNaNs = true
	self.solvers:insert(solver)
	self.exitTime = duration
	self.running = true

	local oldupdate = solver.update
	solver.update = function(self, ...)
		local varname = 'U H'
		local var = assert(self.displayVarForName[varname], "couldn't find "..varname)
		local ymin, ymax, yavg = self:calcDisplayVarRangeAndAvg(var)
	
		local out = table{self.t, ymin, yavg, ymax}:concat'\t'
		print(out)
		f:write(out,'\n')
		
		return oldupdate(self, ...)
	end
end

function App:requestExit()
	App.super.requestExit(self)
os.exit()	-- exit early in Windows to avoid a driver crash
end	

xpcall(function()
	local app = App()
	app:run()

-- work around Windows driver crash
end, function()
	os.exit(1)
end)

f:close()
--unistd.chdir(rundir)

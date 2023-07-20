#!/usr/bin/env luajit
--[[
this will run some order of accuracy tests on different configurations 


--]]

local ffi = require 'ffi'
require 'ffi.c.stdlib'	-- free

local clnumber = require 'cl.obj.number'
local template = require 'template'
local sdl = require 'ffi.sdl'
local class = require 'ext.class'
local table = require 'ext.table'
local range = require 'ext.range'
local fromlua = require 'ext.fromlua'
local tolua = require 'ext.tolua'
local math = require 'ext.math'
local path = require 'ext.path'
local string = require 'ext.string'
local matrix = require 'matrix'
local gnuplot = require 'gnuplot'
local unistd = require 'ffi.c.unistd'

local rundirp = unistd.getcwd(nil, 0)
local rundir = ffi.string(rundirp)
ffi.C.free(rundirp)

-- in case it's needed
local resultsDir = 'results'
path(rundir..'/'..resultsDir):mkdir()

-- from here on the require's expect us to be in the hydro-cl directory
-- I should change this, and prefix all hydro-cl's require()s with 'hydro-cl', so it is require()able from other projects
unistd.chdir'../..'

for k,v in pairs(require 'tests.util') do _G[k] = v end

-- I guess I'll set the HYDROCL_ENV variable to override the sys=console cmdline
-- I could alternatively set the cmdline global ... but for some reason I named this test's cmdline var as 'cmdline' also, so the names collide
--  I could reuse this 'cmdline' as a global to get around that, and set its cmdline.sys='console'
--  or I could rename this 'local cmdline' to 'local testcmdline' and declare a separate global "cmdline={sys='console'}"
-- TODO looks like my own setenv always returns 0, even if it fails to write ...
assert(0 == ffi.C.setenv('HYDROCL_ENV', 'sys="console",'..(os.getenv'HYDROCL_ENV' or ''), 0))

local cmdline = require 'ext.cmdline'(...)

-- which problem to use
local problemName = cmdline.init or 'advect wave'
--local problemName = cmdline.init or 'Sod'

-- don't use cached results <-> regenerate results for selected tests
local nocache = cmdline.nocache

-- for the first configuration, run it at highest resolution, plot it with the exact solution, and quit
local plotCompare = cmdline.compare

-- for the first configuration, plot error vs time from 0 to duration
local plotErrorHistory = cmdline.history

-- exclusive with 'compare': don't use exact, instead use exponential regression
local uselin = cmdline.uselin

local schemeCfgs = require 'ext.fromlua'(path(rundir..'/schemes.lua'):read())

local problems = {}

problems['advect wave'] = {
	configurations = outer(
		{
			{
				eqn = 'euler',
				initCond = 'advect wave',
			}
		},
			-- final error at n=1024 on the right:
		schemeCfgs
	),

	-- TODO make sure solver_t->init_v0x == 1/duration and solver_t->maxs.x - mins.x == 2
	-- otherwise, for durations t=100 and t=1 the results look close enough to the same
	-- or just use what the Mara demo had:
	duration = .1,
}

problems.Sod = {
	-- copy of the above problem ... maybe put somewhere else
	configurations = outer(
		{
			{
				eqn = 'euler',
				initCond = 'Sod',
			}
		},
			-- final error at n=1024 on the right:
			-- (these numbers are all for duration=.1)
		schemeCfgs
	),

	duration = .1,
}

local problem = problems[problemName]

local testdatas = table()
local errorsForConfig = table()
local errorNames = table()

-- for history and compare, this is the size used
--local singleSize = 1024
local singleSize = 64

local dim = 1
local sizes = plotCompare 
	and table{singleSize} 
	or range(3,10):map(function(x) return 2^x end)

if cmdline.time then sizes = table{sizes:last()} end

-- keep C typenames unique
-- but this messes with the cache names as well
-- it would be good if we could dynamically resize the solvers and reuse them ... since size isn't baked into opencl anymore
local allnames = {}

for _,cfg in ipairs(problem.configurations) do
	cfg = table(cfg)

	local destName = string.trim(tolua(cfg):match('^{(.*)}$'):gsub('%s+', ' '):gsub('"', ''))
print(destName)

	local destFilename = destName
		:gsub('/', '')
		:gsub('{', '(')
		:gsub('}', ')')
	
	cfg.dim = dim
	cfg.cfl = .6
	cfg.coord = 'cartesian'

	--[[
	data cached per-test:
	size[size]
		.xs[index] 
		.ys[index] 
	--]]
	local testdata
	local srcfn = rundir..'/'..resultsDir..'/'..destFilename..'.lua'
	local srcfiledata = path(srcfn):read()
	if srcfiledata then
		testdata = fromlua(srcfiledata)
	end
	testdata = testdata or {}
	do	--if nocache or not testdata.size then
		testdata.size = testdata.size or table()
		for _,size in ipairs(sizes) do
			
			testdata.size[size] = testdata.size[size] or table()
			testdata.size[size].errorsForTime = setmetatable(testdata.size[size].errorsForTime or {}, table)
			
			local startTime, endTime
					
			if nocache 
			or not testdata.size[size].xs
			or not testdata.size[size].ys
			or not testdata.size[size].deltaTime
			-- or either no exact or uselin
			then
				cfg = table(cfg)
				cfg.gridSize = {size}
print()
print(size)
print()	

				testdata.size[size].ts = table()

				local duration = tonumber(problem.duration) or error("expected problem.duration")
				
				local function getSolverGraph(solver)
					local err, UBuf = solver:calcExactError(1)

					-- now in solver.reduceBuf
					local numCells = solver.numCells
					local numGhost = solver.numGhost
					-- now in ptr
					
					local xs = range(numCells-2*numGhost):mapi(function(i)
						return (i-.5) * solver.solverPtr.grid_dx.x + solver.solverPtr.mins.x
					end)
					local ys = range(numCells-2*numGhost):mapi(function(i)
						return UBuf[i+numGhost-1].rho
					end)
					local exact
					if not uselin then
						-- TODO this only compares the first value, while 'testAccuracy' cmdline option compares all (integratable variable) state values
						exact = xs:map(function(x)
							return (solver.eqn.initCond:exactSolution(solver.t, x, 0, 0))
						end)
					end
					return xs, ys, exact, assert(err)
				end

				-- some hack trickery: reset the cache names for the opencl prefixes
				-- TODO a better fix would be to just get the live resizing working again
				-- or TODO TODO a better better one would be running the App as a separate process, so the typenames wouldn't differ per cl code, and caching would work fine
				for i=1,math.huge do
					local found
					local function try(name) 
						if allnames[name] then 
							found = true 
							allnames[name] = nil 
						end
					end
					try('initState_'..i)
					try('common_'..i)
					try('solver_'..i)
					try('boundary_'..i)
					if not found then break end
				end
				local App = class(require 'hydro.app')
				App.allnames = allnames
				function App:setup(clArgs)
					cfg.app = self
					local solver = require('hydro.solver.'..cfg.solver)(cfg)
					self.solvers:insert(solver)
					self.exitTime = duration
					self.running = true
					
					local oldupdate = solver.update
					solver.update = function(...)
						local xs, ys, exact, err = getSolverGraph(solver)
						testdata.size[size].ts:insert(solver.t)
						testdata.size[size].errorsForTime:insert(err)
						return oldupdate(...)
					end
					
					startTime = os.clock()
				end
				
				function App:requestExit()
					endTime = os.clock()	-- track time taken

					App.super.requestExit(self)
				
					-- now compare the U buffer to the exact 
					assert(#self.solvers == 1)
					local solver = self.solvers[1]
					local xs, ys, exact, err = getSolverGraph(solver)
					testdata.size[size].xs = xs
					testdata.size[size].ys = ys
					testdata.size[size].exact = exact
					testdata.size[size].error = err
					testdata.size[size].errorsForTime:insert(err)
					testdata.size[size].ts:insert(solver.t)
				end	
				
				local app  = App()
				app:run()
			end
			
			local xs = setmetatable(assert(testdata.size[size].xs), table)
			local ys = setmetatable(assert(testdata.size[size].ys), table)
			local exact = not uselin and assert(testdata.size[size].exact)
			
			if uselin then
				-- just use log/log regression to estimate where the best would be
				-- technically this won't be best, because most our samples tend to flatten at the bottom
				-- and in that case, the regression will point to a less accurate place than where it would if a smaller subset was used 
				local xavg = xs:sum() / #xs
				local yavg = ys:sum() / #ys

				local b1 = range(#xs):map(function(i)
					return (xs[i] - xavg) * (ys[i] - yavg)
				end):sum() / range(#xs):map(function(i)
					return (xs[i] - xavg)^2
				end):sum()
				local b0 = yavg - b1 * xavg
				exact = xs:map(function(x)
					return b0 + b1 * x
				end)
			end

			if plotCompare then -- plotting immediately
				gnuplot{
					output = rundir..'/compare-graphs.svg',
					terminal = 'svg size 2000,1500 background rgb "white"',
					style = 'data lines',
					data = {xs, ys, exact},
					{using = '1:2', title=''..size},
					exact and {using = '1:3', title='exact'} or nil,
				}
				os.exit()
			end

			if path'stop':exists() then path'stop':remove() os.exit(1) end
			if endTime and startTime then
				testdata.size[size].deltaTime = endTime - startTime
			end	
			if cmdline.time then print('time: '..testdata.size[size].deltaTime) end
		end
	end
	testdata.name = destName
	path(srcfn):write(tolua(testdata))
local errors = sizes:map(function(size) return testdata.size[size].error end)
print('error:',table.last(errors))
	errorsForConfig:insert(errors)
	errorNames:insert(destName)
	testdatas:insert(testdata) 
end

if cmdline.time then return end

-- [[ plot errors per grid dx
gnuplot(
	table({
		output = rundir..'/'..resultsDir..'/results '..problemName..'.svg',
		terminal = 'svg size 4000,3000 background rgb "white"',
		style = 'data linespoints',
		log = 'xy',
		xlabel = 'dx',
		ylabel = 'L1 error',
		key = 'left Left reverse samplen 2 spacing .9',
		data = table{
			sizes:map(function(x) return 2/x end),	-- this assumes the domain is -1,1
		}:append(errorsForConfig),
	},
	errorNames:mapi(function(name,i)
		return {using = '1:'..(i+1), title=name}
	end)
))
--]]

-- [[ plot error histories
do
	local graphdir = rundir..'/graphs '..problemName
	path(graphdir):mkdir()
	local data = table.append(testdatas:map(function(testdata)
		return table.append(sizes:map(function(size)
			local sizedata = testdata.size[size]
			local tdata = table(sizedata.ts)
			-- throw away any OOB integration steps
			for i=1,#tdata do
				if tdata[i] > 2*problem.duration then tdata[i] = math.nan end
			end
			return {
				tdata,
				sizedata.errorsForTime,
			}
		end):unpack())
	end):unpack())
	local usings = table.append(testdatas:mapi(function(testdata,i)
		return sizes:mapi(function(size,j)
			local base = 1 + 2 * (j-1 + #sizes * (i-1))
			return {using=base..':'..(base+1), title=testdata.name..', size='..size}
		end)
	end):unpack())
	gnuplot(table({
		output = graphdir..'/error-history.svg',
		terminal = 'svg size 4000,3000 background rgb "white"',
		style = 'data linespoints',
		log = 'xy',
		xlabel = 'time',
		ylabel = 'L1 error',
		nokey = true,	-- can't seem to make the legend small enough to fit them all
		--key = 'left Left reverse samplen 2 spacing .9',
		data = data,
	}, usings))
	
	-- [[ more practical would be to just produce a bunch of collections, like 10 or 20 at a time
	-- using# = 1 + (graphIndex-1) + #sizes * (sizeIndex-1)
	-- compare all graphs of highest resolution
	gnuplot(table({
		output = graphdir..'/error-history-size='..sizes:last()..'.svg',
		terminal = 'svg size 2000,1500 background rgb "white"',
		style = 'data linespoints',
		log = 'xy',
		xlabel = 'time',
		ylabel = 'L1 error',
		key = 'left Left reverse',
		data = data,
	}, range(#sizes,#usings,#sizes):mapi(function(i) return usings[i] end)))
	-- compare each graph's different sizes
	for i=1,#usings,#sizes do
		gnuplot(table({
			output = graphdir..'/error-history-using='..usings[i].title..'.svg',
			terminal = 'svg size 2000,1500 background rgb "white"',
			style = 'data linespoints',
			log = 'xy',
			xlabel = 'time',
			ylabel = 'L1 error',
			key = 'left Left reverse',
			data = data,
		}, usings:sub(i,i+#sizes-1)))
	end
	--]]
end
--]]

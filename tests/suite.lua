#!/usr/bin/env luajit

-- run a suite
-- run this from the test working dir, i.e. hydro/tests/<wherever>/
-- run this as ../suite.lua 
-- TODO instead, maybe run this in the local dir, with the test dir name as the param?  nah, because each test shows/does dif things
--  how about TODO instead, put this in a subdir called "common" or "util" or something?

local ffi = require 'ffi'
local table = require 'ext.table'
local path = require 'ext.path'

-- save the cwd and chdir to ../..
local rundir = path:cwd()

--local resultsDir = 'results'
--path(rundir..'/'..resultsDir):mkdir()
path'../..':cd()

-- set this global to have hydro run in console mode
-- working on doing this cleaner...
cmdline = {
	showfps = true,
	sys = 'console',
	exitTime = 10,
}

local HydroApp = require 'hydro.app':subclass()

-- share across all app instances
HydroApp.allnames = {}

local configs = dofile(rundir..'/config.lua')	-- could use 'require' but there is another in the newly added package.path
for _,config in ipairs(configs) do
	local resultFile = rundir..'/results-'..config.name..'.txt'
	if not path(resultFile):exists() then
		
		-- another TODO for hydro.app ... reuse names for matching ctypes
		local data = table()

		function HydroApp:setup(args)
			print('running config: '..config.name)
			
			local solverClassName = config.solverClassName
			local solverClass = require(solverClassName)
			local args = table(config.solverArgs)
			args.app = self
			local solver = solverClass(args)
			
			local oldUpdate = solver.update
			function solver:update(...)
				oldUpdate(self, ...)

				local row = table()
				row:insert(self.t)
				
				if config.trackVars then
					for _,varName in ipairs(config.trackVars) do
						local var = assert(self.displayVarForName[varName])
						local component = self.displayComponentFlatList[var.component]
						local vectorField = self:isVarTypeAVectorField(component.type)
						local valueMin, valueMax, valueAvg = self:calcDisplayVarRangeAndAvg(var, vectorField and component.magn or nil)
						row:append{valueMin, valueAvg, valueMax}
					end
				end
				
				data:insert(row)
			end

			self.solvers:insert(solver)
		end

		function HydroApp:requestExit()
			path(resultFile):write(
				'#t\t'..table.mapi(config.trackVars, function(varName)
					return varName..' min\t'
						..varName..' avg\t'
						..varName..' max'
				end):concat'\t'..'\n'
				..data:mapi(function(l) 
					return table.concat(l, '\t') 
				end):concat'\n'
				..'\n'
			)
			HydroApp.super.requestExit(self)
		end

		HydroApp():run()
	end
end

path(rundir):cd()

print('configs done -- plotting...')
--os.execute('gnuplot plot.gnuplot')
dofile'../plot.lua'

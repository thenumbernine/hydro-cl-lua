#!/usr/bin/env lua
local tolua = require 'ext.tolua'
local table = require 'ext.table'
local os = require 'ext.os'
local io = require 'ext.io'
local string = require 'ext.string'
local lfs = require 'lfs'

local function updateIndex()
	local fs = table()
	for f in os.listdir'.' do
		if f:match'%.txt$'
		and f ~= 'README.txt'
		then
			fs:insert(f)
		end
	end
	fs:sort()

	local ts = table{[[
<!doctype html>
<html>
	<head>
		<meta charset='utf8'/>
		<title>master</title>
		<style>
table {
    border : 1px solid black;
    border-collapse : collapse;
}
table td, table th {
    border : 1px solid black;
	padding: 4px;
}
		</style>
	</head>
	<body>
		<table>
			<tr>
				<th>name</th>
				<th>time</th>
				<th>avg</th>
				<th>min</th>
				<th>max</th>
				<th>stddev</th>
			</tr>
]]
	}
	for _,f in ipairs(fs) do
		local last = string.split(
			string.split(
				string.trim(io.readfile(f)),
				'\n'
			):last(),
			'%s+'
		)
		-- last = {time, avg, min, max, stddev}
		ts:insert('<tr>'
			..'<td>'
				..io.getfileext(f):gsub('_', ' ')
			..'</td>'
			..last:mapi(function(l) 
				return '<td>'..l..'</td>'
			end):concat()
			..'</tr>')
	end
	ts:insert(
[[
		</table>
	</body>
</html>
]]
	)
	io.writefile('index.html', ts:concat'\n')
end
	

local DIR = lfs.currentdir()

for cfg in coroutine.wrap(function()
	-- TODO vary flux as well?

	for cfgFluxLimiter in coroutine.wrap(function()
		coroutine.yield{
			fluxLimiter = 'donor cell',
			name = '',
		}
	
		-- wow, for the rmin=cylinderRMin boundary, superbee kept it balanced
		-- otherwise, for rmin=none fluxLimiter=superbee it diverges slowly, and for fluxLimiter=donor cell it also diverges slowly
		coroutine.yield{
			fluxLimiter = 'superbee',
			name = '_superbee',
		}
	end) do

-- [[ cartesian
		coroutine.yield{
			coord = 'cartesian',
			mins = '{-1, -1, -1}',
			maxs = '{1, 1, 1}',
			gridSize = '{64, 64, 1}',
			boundary = "{xmin='freeflow', xmax='freeflow', ymin='freeflow', ymax='freeflow', zmin='freeflow', zmax='freeflow'}",
			vectorComponent = 'cartesian',
			fluxLimiter = cfgFluxLimiter.fluxLimiter,
			name = 'cartesian'..cfgFluxLimiter.name,
		}
--]]


		for cfgVectorComponent in coroutine.wrap(function()
	-- [[
			coroutine.yield{
				vectorComponent = 'holonomic',
				name = 'holonomic'..cfgFluxLimiter.name,
			}
			coroutine.yield{
				vectorComponent = 'anholonomic',
				name = 'anholonomic'..cfgFluxLimiter.name,
			}
	--]]		
			coroutine.yield{
				vectorComponent = 'cartesian',
				name = 'cartesian'..cfgFluxLimiter.name,
			}
		end) do
			-- [[ cylinder, rmin == 0
			for cfgBoundary in coroutine.wrap(function()
				-- [=[	boundary rmin == none
				coroutine.yield{
					boundary = "boundary={xmin='none', xmax='freeflow', ymin='periodic', ymax='periodic', zmin='freeflow', zmax='freeflow'}",
					name = cfgVectorComponent.name..'_boundary_none',
				}
				--]=]
				-- [=[	boundary rmin == cylinderRMin
				coroutine.yield{
					boundary = "boundary={xmin='cylinderRMin', xmax='freeflow', ymin='periodic', ymax='periodic', zmin='freeflow', zmax='freeflow'}",
					name = cfgVectorComponent.name..'_boundary_cylinderRMin',
				}
				--]=]	
			end) do
				coroutine.yield(table(cfgVectorComponent, cfgBoundary, {
					coord = 'cylinder',
					mins = '{0, 0, -1}',
					maxs = '{.5, 2*math.pi, 1}',
					gridSize = '{32, 128, 1}',
					fluxLimiter = cfgFluxLimiter.fluxLimiter,
					name = 'cylinder_rmin_eq_0_'..cfgBoundary.name,
				}):setmetatable(nil))
			end
			--]]
		end
	end

--[[ cylinder, rmin != 0
local coord = 'cylinder'
local mins = '{.5, 0, -1}'
local maxs = '{1, 2*math.pi, 1}'
local gridSize = '{32, 128, 1}'
local boundary = "boundary={xmin='freeflow', xmax='freeflow', ymin='periodic', ymax='periodic', zmin='freeflow', zmax='freeflow'}"
local name = 'cylinder_rmin_gt_0'
--]]

-- TODO sphere
-- TODO sphere_sinh_radial

end) do

	assert(lfs.chdir'../..')

	-- forward args
	local cmd = './run.lua sys=console "config='..DIR..'/config.lua" exitTime=1'
	cmd = cmd .. ' "trackvars=U Pi"'
	cmd = cmd .. ' coord='..cfg.coord
	cmd = cmd .. ' vectorComponent='..cfg.vectorComponent
	cmd = cmd .. ' "mins='..cfg.mins..'"'
	cmd = cmd .. ' "maxs='..cfg.maxs..'"'
	cmd = cmd .. ' "gridSize='..cfg.gridSize..'"'
	cmd = cmd .. ' "fluxLimiter='..cfg.fluxLimiter..'"'
	cmd = cmd .. ' "'..cfg.boundary..'"'

	-- what's the difference between plotOnExit and plot1DOnExit?
	cmd = cmd .. ' "plotOnExit='..DIR..'/'..cfg.name..'.svg"'
	-- columns are t, avg, min, max, stddev
	cmd = cmd .. ' "plotOnExit_savedata='..DIR..'/'..cfg.name..'.txt"'

	print()
	print(cmd)
	print(os.execute(cmd))

	assert(lfs.chdir(DIR))
	
	-- update our global matrix and write out
	updateIndex()
end



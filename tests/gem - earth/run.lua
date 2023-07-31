#!/usr/bin/env luajit
require 'ext'
local exec = require 'exec'	-- TODO put this somewhere everyone can get it
local ffi = require 'ffi'

local cwd = path:cwd()
print('cwd = '..cwd)
path'../..':cd()

for _,coord in ipairs{
	'cartesian',
	'cylinder',
	'sphere',
} do
	for _,vectorComponent in ipairs{
		--'anholonomic',
		'cartesian',
	} do
		-- adjust grid size to produce less rows for coord charts like sphere that have small dx's near the center
		local gridSize = {32,32,32}
		if coord == 'sphere' then gridSize = {16,16,16} end

		local cfl = .1

		-- TODO also vary vectorComponent, gridSize, useCTU?
		local dest = table{
			'out',

		-- [[ matches below
			'coord='..coord,
			'vectorComponent='..vectorComponent,
			'gridSize={'..table.concat(gridSize,',')..'}',
			'cfl='..cfl,
		--]]
		}:concat' '..'.txt'

		-- for how long should this run?
		-- TODO pick coord ... and pick vectorComponent
		local cmd = table{
			'luajit',
			'run.lua',
			'sys=console',
			'exitTime=.5',
			'"trackvars=U E_g mag"',
			'"config='..cwd..'/config.lua"',

		-- [[ matches above
			'coord='..coord,
			'vectorComponent='..vectorComponent,
			'"gridSize={'..table.concat(gridSize,',')..'}"',	-- need ""'s to wrap {}'s ?
			'cfl='..cfl,
		--]]

			'> "'..cwd..'/'..dest..'"',
		}:concat' '
		exec(cmd)
	end
end

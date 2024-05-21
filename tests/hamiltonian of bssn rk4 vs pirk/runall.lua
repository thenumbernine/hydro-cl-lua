#!/usr/bin/env luajit
for dim=1,3 do
	for _,solvername in ipairs{
		'bssnok-fd',
		'bssnok-pirk',
	} do
		for _,initState in ipairs{
			'Minkowski',
			'pure gauge wave',
			--'scalar field',
		} do
			local cmd = 'luajit run.lua '
				..dim..' '
				..('%q'):format(solvername)..' '
				..('%q'):format(initState)
			print(require 'ext.os'.exec(cmd))
		end
	end
end

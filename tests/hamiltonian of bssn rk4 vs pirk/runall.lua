for dim=2,2 do
	for _,solvername in ipairs{
		'bssnok-fd',
		--'bssnok-pirk',
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
			print('>'..cmd)
			print('exec results',os.execute(cmd))
		end
	end
end

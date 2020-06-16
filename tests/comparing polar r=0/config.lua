return dofile'tests/comparing spiral axisymmetric/config.lua':mapi(function(config)
	-- TODO maybe merge these two, since mesh solver doesn't use args.mins 
	config.solverArgs.mins = {0,0,0}
	if config.solverArgs.mesh then
		config.solverArgs.mesh.mins = {0,0,0}
	end
	return config
end)

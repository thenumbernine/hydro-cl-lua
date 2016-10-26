return function(name, vars)
	return 'typedef struct { real ' .. table.concat(vars, ', ') .. '; } '..name..';'
end

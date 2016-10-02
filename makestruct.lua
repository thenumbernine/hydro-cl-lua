return function(name, vars)
	return 'typedef struct { real ' .. vars:concat', ' .. '; } '..name..';'
end

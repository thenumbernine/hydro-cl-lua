local table = require 'ext.table'
--[[
idk where to put this stuff
--]]
return function(env)
	env = env or {}
	local xNames = table{'x', 'y', 'z'}
	env.xNames = xNames

	env.symNames = table{'xx', 'xy', 'xz', 'yy', 'yz', 'zz'}

	local from3x3to6_table = {{1, 2, 3}, {2, 4, 5}, {3, 5, 6},}
	env.from3x3to6 = function(i,j) return from3x3to6_table[i][j] end

	local from6to3x3_table = {{1,1},{1,2},{1,3},{2,2},{2,3},{3,3}}
	env.from6to3x3 = function(i) return table.unpack(from6to3x3_table[i]) end

	env.sym = function(a,b)
		assert(a >= 1 and a <= 3, "tried to index sym with "..tostring(a)..", "..tostring(b))
		assert(b >= 1 and b <= 3, "tried to index sym with "..tostring(a)..", "..tostring(b))
		if a > b then a,b = b,a end
		return xNames[a]..xNames[b]
	end

	return env
end

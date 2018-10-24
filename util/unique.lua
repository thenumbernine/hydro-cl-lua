-- static so that no names overlap across all equations
local allnames = {}
function uniqueName(name)
	--[[ I'm using the base name in my typedefs per-source file
	if not allnames[name] then
		allnames[name] = true
		return name
	end
	--]]
	for i=1,math.huge do
		local try = name..'_'..i
		if not allnames[try] then
			allnames[try] = true
			return try
		end
	end
end
return uniqueName

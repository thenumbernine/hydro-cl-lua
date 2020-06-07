local function roundup(a, b)
	local mod = a % b
	if mod ~= 0 then
		a = a - mod + b
	end
	return a
end

return roundup

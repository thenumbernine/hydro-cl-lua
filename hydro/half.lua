local ffi = require 'ffi'
local half = require 'cl.obj.half'

local function fromreal(x)
	if ffi.sizeof'real' == 2 then
		--ffi.istype('half', x)
		return half.from(x)
	end
	return x
end	
	
local function toreal(x)
	if ffi.sizeof'real' == 2 then
		return half.to(x)
	end
	return x
end

--[[ how well does this work?
--[=[
local function test(x)
	local h = half.to(x)
	local y = half.from(h)
	print(math.abs(y-x), x, ('%x'):format(h.i), y)
end
for x=-10,10 do test(x + 1e-5) end
for x=-10,10 do test(x - 1e-5) end
for x=-10,10 do test(10^x) end
for x=-10,10 do test(-10^x) end
test(1/0)
test(-1/0)
test(0/0)
test(-0/0)
--]=]

local math = require 'ext.math'
-- numbers that don't map back and forth: 16-bit nans and 0
for i=0,65536 do
	local h = ffi.new'half'
	h.i = i
	local x = half.from(h)
	assert(type(x) == 'number')
	b32.f = x
	local h2 = half.to(x)
	if i ~= h2.i then
		print(('%x\t%x\t%x'):format(i, h2.i, b32.i), x)
		assert(not math.isfinite(x) or h.i == h2.i)
	end
end
os.exit()
--]]

return {
	from = half.from,
	to = half.to,
	
	fromreal = fromreal,
	toreal = toreal,
}

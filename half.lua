local bit = require 'bit'
local ffi = require 'ffi'
local math = require 'ext.math'

ffi.cdef[[
typedef union {
	unsigned int i;
	float f;
	unsigned char ptr[4];
	struct {
		unsigned int mant : 23;
		unsigned int exp : 8;
		unsigned int sign : 1;
	};
} float32bits_t;
]]
assert(ffi.sizeof'float32bits_t' == 4)

local b32 = ffi.new'float32bits_t'
local function tohalf(x)	-- converts to 'half' defined in cl/obj/env.lua
	b32.f = x

	local h = ffi.new'half'
	h.sign = 0
	if b32.exp == 0xff then		-- nan / inf
		h.mant = bit.rshift(b32.mant, 13)
		if b32.mant ~= 0 then
			h.mant = bit.bor(h.mant, 0x200)
		end
		h.exp = 0x1f		-- nan / inf
	else
		local exp = b32.exp - 0x70
		if exp >= 0x1f then			-- exp overflow -> +-inf
			h.mant = 0
			h.exp = 0x1f
		elseif exp < -10 then	-- zero
			h.mant = 0
			h.exp = 0
		else
			if exp <= 0 then		-- exp in [-10,0]
				local negexp = 13 - exp	-- negexp in [13,23]
				local c = bit.bor(b32.mant, 0x800000)
				h.i = bit.rshift(c, negexp+1)
				if bit.band(c, bit.lshift(1, negexp)) ~= 0 
				and bit.band(c, (bit.lshift(3, negexp) - 1)) ~= 0 
				then
					h.i = h.i + 1
				end
			else
				h.mant = bit.rshift(b32.mant, 13)
				h.exp = exp
				if bit.band(b32.mant, 0x1000) ~= 0 
				and bit.band(b32.mant, 0x2fff) ~= 0 
				then
					h.i = h.i + 1	-- add through mant & exp
				end
			end
		end
	end
	h.sign = b32.sign
	return h
end

local function fromhalf(x)
	local exp = bit.rshift(bit.band(x.i, 0x7c00), 10)
	local mant = bit.lshift(bit.band(x.i, 0x3ff), 13)

	if exp == 0x1f then
		exp = 0xff
		if mant ~= 0 then
			mant = bit.bor(mant, 0x400000)
		end
	elseif exp == 0 then
		if mant ~= 0 then
			exp = exp + 1
			while bit.band(mant, 0x7f800000) == 0 do
				mant = bit.lshift(mant, 1)
				exp = exp - 1
			end
			mant = bit.band(mant, 0x7fffff)
			exp = exp + 0x70
		end
	else
		exp = exp + 0x70
	end

	b32.mant = mant
	b32.exp = exp
	b32.sign = x.sign
	
	--[[
	looks like tonumber(float32) doesn't work for all possible float binary representations ... like 0xfff80000
	Lua print() and Lua tonumber() treat values as null, but printf() treats them like floats
	
	7f800000 = 0 : 111 1111 1000 0 : 000 0000 0000 0000 0000 = inf
	ff800000 = 1 : 111 1111 1000 0 : 000 0000 0000 0000 0000 = -inf
	ffc00000 = 1 : 111 1111 1100 0 : 000 0000 0000 0000 0000 = nan
	
	value doesn't convert with tonumber(float32): 
	
	fff80000 = 1 : 111 1111 1111 1 : 000 0000 0000 0000 0000
	--]]
	local y = tonumber(b32.f)
	if not y then y = math.nan end
	--assert(y, "failed to convert for "..('%x'):format(b32.i)..' '..tostring(b32.f))
	return y
end


local function fromreal(x)
	if ffi.sizeof'real' == 2 then
		--ffi.istype('half', x)
		return fromhalf(x)
	end
	return x
end	
	
local function toreal(x)
	if ffi.sizeof'real' == 2 then
		return tohalf(x)
	end
	return x
end

--[[ how well does this work?
--[=[
local function test(x)
	local h = tohalf(x)
	local y = fromhalf(h)
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
	local x = fromhalf(h)
	assert(type(x) == 'number')
	b32.f = x
	local h2 = tohalf(x)
	if i ~= h2.i then
		print(('%x\t%x\t%x'):format(i, h2.i, b32.i), x)
		assert(not math.isfinite(x) or h.i == h2.i)
	end
end
os.exit()
--]]

return {
	from = fromhalf,
	to = tohalf,
	
	fromreal = fromreal,
	toreal = toreal,
}

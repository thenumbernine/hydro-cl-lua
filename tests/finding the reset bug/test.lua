#!/usr/bin/env luajit
require 'ext'
local Image = require 'image'
local matrix = require 'matrix'

local bufnames = {
	'deltaUEigBuf',
	'eigenBuf',
	'fluxBuf',
	'reduceBuf',
	'rEigBuf',
	'UBuf',
	'waveBuf',
}
local imgs = table()
for i=1,2 do
	for j=1,2 do
		for _,name in ipairs(bufnames) do
			local prefix = 'r'..i..'f'..j..'_'..name
			local filename = prefix..'.fits'
			imgs[prefix] = Image(filename)
		end
	end
end

local function unravel(index,size)
	local result = matrix()
	for i,s in ipairs(size) do
		local j = index % s
		index = (index - j) / s
		table.insert(result, j)
	end
	table.insert(result, index)
	return result
end

local function compare(fa,fb)
	local A = imgs[fa]
	local B = imgs[fb]
	
	local imgsize = matrix{A:size()}
	assert(imgsize == matrix{B:size()})
	local size = matrix{36,36,select(2, imgsize:unpack())}
	
	local chs = table()
	local errs = table()
	local diffs 
	print('comparing '..fa..' with '..fb..'...')
	for i=0,size:prod()-1 do
		local x1, x2 = A.buffer[i], B.buffer[i]
		if x1 ~= x2 then 
			diffs = true
			local ch = select(4,unravel(i,size):unpack())
			--print(unravel(i,size), x1, x2) 
			chs[ch] = (chs[ch] or 0) + 1
			errs[ch] = (errs[ch] or 0) + math.abs(x1 - x2)
		
		end 
	end
	if diffs then
		print('...differs!')
		print(tolua(chs, {indent=true}))
		print(tolua(errs, {indent=true}))
	end
end

for i=1,2 do
	for _,name in ipairs(bufnames) do
		local n1 = 'r1f'..i..'_'..name
		local n2 = 'r2f'..i..'_'..name
		compare(n1, n2)
	end
end

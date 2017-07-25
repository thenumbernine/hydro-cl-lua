#!/usr/bin/env luajit
require 'ext'
local Image = require 'image'
local matrix = require 'matrix'

local imgs = table{
	r1f1_U = Image'r1f1_U.fits',
	r1f2_U = Image'r1f2_U.fits',
	r2f1_U = Image'r2f1_U.fits',
	r2f2_U = Image'r2f2_U.fits',
	r1f1_wave = Image'r1f1_wave.fits',
	r1f2_wave = Image'r1f2_wave.fits',
	r2f1_wave = Image'r2f1_wave.fits',
	r2f2_wave = Image'r2f2_wave.fits',
}

local imgsize = matrix{imgs.r1f1_U:size()}
local size = matrix{36,36,select(2, imgsize:unpack())}

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

compare('r1f1_U', 'r2f1_U')
compare('r1f1_wave', 'r2f1_wave')
compare('r1f2_U', 'r2f2_U')
compare('r1f2_wave', 'r2f2_wave')

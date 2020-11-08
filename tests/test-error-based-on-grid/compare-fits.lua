#!/usr/bin/env luajit
require 'ext'
local Image = require 'image'
local imgs =
table.wrapfor(os.listdir'.'):mapi(function(x)
	return x[1]
end)	-- filenames
:filter(function(fn)
	return fn:match'%.fits$'
end)	-- fits filenames
:sort(function(a,b)
	local pat = 'U%-(%d+)%.fits'
	return tonumber(a:match(pat)) < tonumber(b:match(pat))
end)	-- filenames sorted smallest to largest
:mapi(function(fn)
	return Image(fn)
end)
local function calcValue(im, x, y)
	-- [[ rho
	return im.buffer[0 + im.channels * (x + im.width * y)]
	--]]
end

local errors = table()
local volumes = table()

local numghost = 2
for i=1,#imgs-1 do
	local lo = imgs[i]
	local hi = imgs[i+1]
	local lw, lh, hw, hh = lo.width-2*numghost, lo.height-2*numghost, hi.width-2*numghost, hi.height-2*numghost
	assert(lw * 2 == hw)
	assert(lh * 2 == hh)
	local diffnorm = 0
	for y=0,hw-1 do
		for x=0,hh-1 do
			local vhi = calcValue(hi,x,y)
			local vlo = calcValue(lo, bit.rshift(x,1)+numghost, bit.rshift(y,1)+numghost)
			local diff = vhi - vlo
			diffnorm = diffnorm + diff * diff
		end
	end
	diffnorm = math.sqrt(diffnorm) / (hw * hh)
	local volume = hw * hh
	volumes:insert(volume)
	errors:insert(diffnorm)
	print(volume, diffnorm)
end

print'plotting...'
local gnuplot = require 'gnuplot'
gnuplot{
	output = 'convergence.png',
	data = {
		volumes,
		errors,
	},
	style = 'data linespoints',
	log = 'xy',
	xtitle = 'volume',
	ytitle = 'error',
	{using = '1:2', title = 'error vs volume'}
}

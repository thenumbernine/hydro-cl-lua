#!/usr/bin/env luajit
require 'ext'
local matrix = require 'matrix'
local _ = matrix.index
local ms = table()
local fns = table()
--[[ out of order
for fn in file:dir() do
--]]
-- [[ in order
for _,fn in ipairs{
	'euler 2d gmres err - 16x16.txt',
	'euler 2d gmres err - 32x32.txt',
	'euler 2d gmres err - 64x64.txt',
	'euler 2d gmres err - 128x128.txt',
	'euler 2d gmres err - 256x256.txt',
	'euler 3d gmres err - 16x16x16.txt',
	'euler 3d gmres err - 32x32x32.txt',
} do
--]]
	if select(2, io.getfileext(fn)) == 'txt' then
		fns:insert(fn)
		local ls = file(fn):read():trim():split'\n'
		local firstline = ls:remove(1)
		assert(firstline:sub(1,1) == '#', "expected first line to be a comment, found:\n"..firstline)
		local m = matrix(
			ls:map(function(l)
				return l:split'%s+':map(function(w)
					return (assert(tonumber(w), "tonumber failed for "..w))
				end)
			end)
		):transpose()
		local iters = m[3]
		for line,iter in ipairs(iters) do 
			assert(iter == 1, "found non-one iteration in "..fn..' line '..line..': '..iters[i])
		end
		ms:insert(m[1])
		ms:insert(m[2])
	end
end

local sigma = 100
local kexp = matrix(range(-30, 30)):map(function(x) return math.exp(-(x/sigma)^2) end)

function matrix.kernel(m, k, normalize)
	local normalization = normalize and k:sum() or 1
	local r = matrix(range(#k))
	return (m:size() - #k):lambda(function(i)
		return m(r+i-1) * k
	end)
end

ms = ms:map(function(m)
	return m:kernel(kexp, true)
end)

require 'gnuplot'(table({
	output = 'results.png',
	style = 'data lines',
	log = 'y',
	format = {y = '%.3e'},
	data = ms,
}, 
	fns:map(function(fn,i)
		return {using=(2*i-1)..':'..(2*i), title=fn}
	end)
))

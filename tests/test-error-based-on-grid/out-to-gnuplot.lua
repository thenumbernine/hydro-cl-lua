#!/usr/bin/env lua
require 'ext'
local lfs = require 'lfs'
local DIR = lfs.currentdir()
-- transform 'out.txt' into a gnuplot-based file
local outdst = DIR..'/out.txt'
local datafn = DIR..'/data.txt'
local ks
local data = file[outdst]:split'\n':mapi(function(l,i,t)
	-- first line always starts with t=
	if l:match'^t=' then
		local firstcol
		if not ks then
			ks = table()
			firstcol = true
		end	
		local ws = l:split'\t':mapi(function(w)
			local k, v = w:match'([^=]*)=(.*)'
			k = k:gsub('%s', '_')
			local vs = v:match'^%[(.*)%]$'
			if vs then
				k = table{k..'_min', k..'_avg', k..'_max'}:concat'\t'
				v = vs:gsub(' ','\t')
			end
			if firstcol then
				ks:insert(k)
			end
			return v
		end)
		return ws:concat'\t', #t+1
	end
end):concat'\n'..'\n'
file[datafn] = '#'..ks:concat'\t'..'\n'..data

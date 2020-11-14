#!/usr/bin/env lua
require 'ext'
local exec = require 'exec'
local gnuplot = require 'gnuplot'

-- transform 'out.txt' into a gnuplot-based file
local outdst = 'out.txt'
local datafn = 'data.txt'
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

local datafn = 'data.txt'
local cols = file[datafn]:split'\n'[1]:sub(2):split'\t'
for i=2,#cols do
	local col = cols[i]
	gnuplot{
		output = 'out-'..col..'.png',
		style = 'data lines',
		{datafile = datafn, using='1:'..i},
	}
end
for i=2,#cols,3 do
	local col = cols[i]
	local prefix = assert(col:match'(.*)_min$')
	gnuplot{
		savecmds = 'gnuplot-cmds.txt',
		output = 'out-'..prefix..'.png',
		style = 'data lines',
		{datafile = datafn, using='1:'..i..':'..(i+2), title='range', with='filledcurve', fillstyle='transparent solid 0.5', linetype='rgb "#ff0000"'},
		{datafile = datafn, using='1:'..i, title='min', linetype='rgb "#ff0000"'},
		{datafile = datafn, using='1:'..(i+1), title='avg', linetype='rgb "#ff0000"'},
		{datafile = datafn, using='1:'..(i+2), title='max', linetype='rgb "#ff0000"'},
	}
end

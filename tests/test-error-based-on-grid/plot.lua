#!/usr/bin/env lua
require 'ext'
local exec = require 'exec'
local gnuplot = require 'gnuplot'

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

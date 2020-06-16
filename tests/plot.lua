#!/usr/bin/env lua
require 'ext'
local gnuplot = require 'gnuplot'
for f in io.dir'.' do
	local ident = f:match'^results%-(.*)%.txt$'
	if ident then
		gnuplot{
			terminal = 'png size 2400,1400',
			output = 'results-'..ident..'.png',
			style = 'data lines',
			xlabel = 'time',
			ylabel = 'velocity',
			xrange = {0,10},
			yrange = {0, .55},	-- |v| starts at .5 and decreases
			key = 'left Left reverse',
			{datafile='results-'..ident..'.txt', using='1:2', title='v min'},
			{datafile='results-'..ident..'.txt', using='1:3', title='v avg'},
			{datafile='results-'..ident..'.txt', using='1:4', title='v max'},
		}

		gnuplot{
			terminal = 'png size 2400,1400',
			output = 'diff'..ident..'.png',
			style = 'data lines',
			xlabel = 'time',
			ylabel = 'velocity',
			xrange = {0,10},
			key = 'left Left reverse',
			{datafile='results-'..ident..'.txt', using='1:($2-$3)', title='v min-avg'},
			{datafile='results-'..ident..'.txt', using='1:($4-$3)', title='v max-avg'},
		}
	end
end

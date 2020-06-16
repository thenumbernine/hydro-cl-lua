#!/usr/bin/env lua
require 'ext'
local gnuplot = require 'gnuplot'
for f in io.dir'.' do
	local ident = f:match'^results%-(.*)%.txt$'
	if ident then
		local names = file[f]:split'\n'[1]:sub(2):split'\t'

		local args, argsDiff
		if io.fileexists'gnuplot-config.lua' then
			args, argsDiff = table.unpack(dofile'gnuplot-config.lua')
		end

		gnuplot(table(args, {
			terminal = 'png size 2400,1400',
			output = 'results-'..ident..'.png',
			style = 'data lines',
			xlabel = names[1],
			ylabel = names[2]:sub(1,-5),
			xrange = {0,10},
			key = 'left Left reverse',
			{datafile=f, using='1:2', title=names[2]},
			{datafile=f, using='1:3', title=names[3]},
			{datafile=f, using='1:4', title=names[4]},
		}))

		gnuplot(table(argsDiff, {
			terminal = 'png size 2400,1400',
			output = 'diff'..ident..'.png',
			style = 'data lines',
			xlabel = names[1],
			ylabel = names[2]:sub(1,-5),
			xrange = {0,10},
			key = 'left Left reverse',
			{datafile=f, using='1:($2-$3)', title=names[2]..' - '..names[3]},
			{datafile=f, using='1:($4-$3)', title=names[4]..' - '..names[3]},
		}))
	end
end

#!/usr/bin/env luajit
require 'ext'
local gnuplot = require 'gnuplot'

-- hmm, how to union all data?  map key = time, union rows, insert '-''s where there are no matching records?
local datafns = table()

-- parse trackvars out of output
for fn in file:dir() do
	local base = fn:match'^out (.*)%.txt'
	if base then
		print('processing '..fn)
		local data = table()
		for l in io.lines(fn) do
			local t, min, avg, max = l:match'^t=(%S+)\tU E_g mag=%[(%S+) (%S+) (%S+)%]$'
			if t then
				data:insert(table{t, min, avg, max}:mapi(function(x) return tonumber(x) or '-' end))
			end
		end
		print('got '..#data..' rows')
		if #data > 0 then
			local datafn = 'plotdata '..base..'.txt'
			file(datafn):write(data:mapi(function(row)
				return row:concat'\t'
			end):concat'\n')
			datafns:insert(datafn)
			-- TODO only regen upon request? or nah?
			gnuplot{
				terminal = 'svg size 1024,768',
				output = fn:gsub('%.txt$', '.svg'),
				style = 'data lines',
				xlabel = 't',
				ylabel = '|gravity| (m/s^2)',
				--data = require 'matrix'(data):T(),
				{datafile=datafn, using = '1:2', title='min'},
				{datafile=datafn, using = '1:3', title='avg'},
				{datafile=datafn, using = '1:4', title='max'},
			}
		end
	end
end

do return end
if #datafns > 0 then
	local args = table(
		{
			terminal = 'svg size 1024,768',
			output = 'plot-all.svg',
			style = 'data lines',
			xlabel = 't',
			ylabel = '|gravity| (m/s^2)',
		}, 
		table():append(
			datafns:mapi(function(fn)
				return {
					{datafile=fn, using = '1:2', title=fn..' min'},
					{datafile=fn, using = '1:3', title=fn..' avg'},
					{datafile=fn, using = '1:4', title=fn..' max'},
				}
			end):unpack()
		)
	):setmetatable(nil)
	print(tolua(args))
	gnuplot(args)
end

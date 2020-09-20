#!/usr/bin/env luajit
require 'ext'
local fs = table()
for f in os.listdir'.' do
	if f:match'%.txt' then fs:insert(f) end
end

local maxrows = 200
os.execute('gnuplot -p -e "'..table{
	"set terminal png size 1024,768",
	"set style data lines",
	"set log y",
	"set output 'residual.png'",
	"plot [0:"..maxrows.."] "..fs:map(function(f)
		return "'"..f.."' using 1:2"
	end):concat', ',
	"set output 'x_norm.png'",
	"plot [0:"..maxrows.."] "..fs:map(function(f)
		return "'"..f.."' using 1:3"
	end):concat', ',
	"unset log y",
	"set output 'x_min_max.png'",
	"plot [0:100] "..fs:map(function(f)
		return "'"..f.."' using 1:4:5 with filledcurves fs transparent solid .5"
	end):concat', ',
}:concat'; ')

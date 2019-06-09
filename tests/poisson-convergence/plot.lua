#!/usr/bin/env luajit
require 'ext'
local fs = table()
for f in io.dir'.' do
	if f:match'%.txt' then fs:insert(f) end
end
os.execute('gnuplot -p -e "'..table{
	"set terminal png size 1024,768",
	"set output 'out.png'",
	"set style data lines",
	"set log y",
	"plot [0:100] "..fs:map(function(f)
		return "'"..f.."' using 1:2"
	end):concat', '
}:concat'; ')

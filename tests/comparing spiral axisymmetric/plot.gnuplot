#!/usr/bin/env gnuplot

set terminal png size 2400,1400

set style data lines
set xlabel "time"
set ylabel "velocity"
set key left Left reverse

set output "results.png"
plot "results.txt" using 1:2 title "v min", \
	"results.txt" using 1:3 title "v avg", \
	"results.txt" using 1:4 title "v max"

set output "diff.png"
plot "results.txt" using 1:($2-$3) title "v min-avg", \
	"results.txt" using 1:($4-$3) title "v max-avg"

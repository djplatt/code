set term postscript eps enhanced color
set output "interval_circ.eps"
#set xrange [-6:4]
set grid
plot "interval_circ.txt" with points


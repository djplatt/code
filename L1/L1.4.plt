set term postscript eps enhanced color
set output "L1.4.eps"
set grid
set nokey
set xlabel "q=modulus"
set ylabel "max |L(1,{/Symbol c})|-1/2 log (q) +0.02012"
set xrange [0:2000000]
#
plot "L1.4_max.out" using 1:($2-log($1)/2+0.02012) pt 1 ps .1

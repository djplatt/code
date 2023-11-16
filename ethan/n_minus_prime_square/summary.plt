set terminal pdf
set output summary.pdf
set nokey
set xlabel "$\log_2 N$"
set ylabel "$A(N)$"
set grid
plot "summary.lst" with lines
#,10498/2*(1 + erf(0.225*(x - 19)))
quit

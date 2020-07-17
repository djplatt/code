set terminal postscript enhanced color eps
set output "kevin1.eps"
set grid
plot "kevin_10_1e-4_4.out" title "Q=10" with lines,"kevin_100_1e-4_4.out" title "Q=100" with lines,"kevin_1000_1e-4_4.out" title "Q=1000" with lines,"kevin_10000_1e-4_4.out" title "Q=10000" with lines,"kevin_100000_1e-4_4.out" title "Q=100000" with lines
quit

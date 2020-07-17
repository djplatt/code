set term postscript eps mono
binwidth=0.01
set boxwidth binwidth
bin(x,width)=width*floor(x/width)+binwidth/2.0
set xlabel "Normalised Gap"
set ylabel "Relative Frequency"
set xrange [0:3.5]
set output "zero_gaps.eps"
plot "gaps_30607946000.dat" using (bin($1,binwidth)):(1.0)/74554.58 smooth freq title "near 30,608,996,000","gaps_99146000.dat" using (bin($1,binwidth)):(1.0)/55430.42 smooth freq title "near    100,196,000", "gaps_10946000.dat" using (bin($1,binwidth)):(1.0)/48332.08 smooth freq title "near     11,996,000","gaps_446000.dat" using (bin($1,binwidth)):(1.0)/41049.94 smooth freq title "near      1,496,000"

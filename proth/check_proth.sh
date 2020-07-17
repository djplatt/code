#!/bin/bash
# run check_proth
#
# use 6 cores
#
ncores=6
c=0
# do 2.5e10 h each time.
step=25000000000
# we need h>=10^27/2^45 < 2.85e13
lim=28500000000000
# eliminate those divisible by primes < 14 
for((h0=0;h0<$lim;h0+=$step))
do
   h1=`expr $h0 + $step`
   ./check_proth $HOME/data/proth/proth_${h0}_${h1}.dat > $HOME/data/proth/check_${h0}_${h1}.log &
   c=`expr $c + 1`
   if [ $c -eq $ncores ]
   then
      c=0
      wait
   fi
done
wait

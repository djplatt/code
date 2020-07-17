#!/bin/bash
#
n=52
h0=4500000000000000
h1=`expr $h0 + 4000000000`
time ./proth1.4 $n $h0 $h1 16000 proth_${n}_${h0}_${h1}.dat
time ./check_proth proth_${n}_${h0}_${h1}.dat
# rm proth_$n_$h0_$h1.dat

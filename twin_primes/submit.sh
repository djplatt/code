#!/bin/bash
st=0
its=32768
bigstep=5242880000000000
for((n=0;n<24;n++))
do
   qsub -joe -v "st=${st},its=${its}" do_twin_primes.sh
   st=$(($st + $bigstep))
done

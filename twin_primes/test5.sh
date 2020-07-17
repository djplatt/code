#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
cd $HOME/code/twin_primes
rm -f results5.dat
st=100000000000000000
cs=16
gap=10000000000
its=16
biggap=$(($its * $gap))
for((n=0;n<cs;n++))
do
   en=$(($st + $biggap))   
   ./twin_primes $st $en $gap >> results5.dat &
   st=$en
done
wait
date

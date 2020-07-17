#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
cd $HOME/code/twin_primes
rm -f results2.dat
st=100000000000000000
for((m=0;m<64;m++))
do
for((n=0;n<16;n++))
do
en=$(($st + 10000000000))
echo "$st $en $($HOME/primesieve-5.4.1/primesieve -c2 -t1 $st $en | awk '/Twin/ {print $4}')" >> results2.dat &
st=$en
done
wait
done
date

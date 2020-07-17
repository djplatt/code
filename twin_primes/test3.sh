#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
cd $HOME/code/twin_primes
rm -f results3.dat
st=100000000000000000
for((n=0;n<1024;n++))
do
en=$(($st + 10000000000))
echo "$st $en $($HOME/primesieve-5.4.1/primesieve -c2 -t16 $st $en | awk '/Twin/ {print $4}')" >> results3.dat
st=$en
done
date

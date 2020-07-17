#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
ncores=16
n=0
nn=0
for((st=130646000;;st+=2100000))  
do
    $HOME/code/chris/mpfi_test 200 $HOME/data/first_billion/zeros_${st}.dat 72795 100 1.64e19 > $HOME/data/chris/mpfi_test_${st}.out &
    n=`expr $n + 1`
    nn=`expr $nn + 1`
    if [ $n -eq $ncores ]
    then
       n=0
       wait
    fi
    if [ $nn -eq 256 ]
    then
       wait
       break
    fi
done
date

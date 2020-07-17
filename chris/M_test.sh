#!/bin/bash
#
#PBS -l walltime=11:45:00,nodes=1:ppn=16
#PBS -q veryshort
date
hostname
ncores=16
mkdir -p $HOME/data/chris/all_asearch_${alpha}/$om
n=0
for f in $(ls $HOME/data/zeros/zeros_*.dat)
do
    ${HOME}/code/chris/M_sum 200 ${f} ${om} 10000000000 $alpha >> ${HOME}/data/chris/all_asearch_${alpha}/${om}/${n}.out &
    n=`expr $n + 1`
    if [ $n -ge $ncores ]
    then
        n=0
        wait
    fi
done
wait
grep "test on" ${HOME}/data/chris/all_asearch_${alpha}/$om/*.out | awk '//{print $19 " " $20 " " $21}' | $HOME/code/chris/arb_sum
date

#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
for om in 7279513355
   mkdir $HOME/data/chris/$om
   for f in $(ls $HOME/data/first_billion/zeros_*.dat)
   do
       $HOME/code/chris/M_sum 200 $f $om 10000000 1.64e19 > $HOME/data/chris/$om/${f}.out &
   done
   wait
done
date

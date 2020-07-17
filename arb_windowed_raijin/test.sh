#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
#
#---------------------------------------------
# you would edit this section
# 
#
#---------------------------------------------
#   everything from here on is standard
#
date
if [ -z "$step" ]
   then
   echo "must specify step"
   exit
fi
hostname
nthreads=16
bigstep=`expr $step \* 20`
start=49000000000
for((c=0;c<$nthreads;c++))
do
   $HOME/code/arb_windowed_raijin/arb_zeta 200 $start 20 $step > $HOME/code/arb_windowed_raijin/results/out_16_${c}.log &
   start=`expr $start + $bigstep`
done
wait
date

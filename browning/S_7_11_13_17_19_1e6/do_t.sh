#!/bin/bash
#
#PBS -l walltime=240:00:00,nodes=1:ppn=16
#PBS -q medium
date
hostname
cd $HOME/code/browning/S_7_11_13_17_19_1e6/
for t in $(cat $dfile)
do
   ../do_t $t 1000000 > do_t_${t}.log &
done
wait
date

#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
ncores=16
zd="${HOME}/data/first_billion"
n=0
for((core=0;core<${ncores};core++))
do
   om0=$((138000 + $core * 37500))
   om1=$(($om0 + 37500))
   ${HOME}/code/chris/big_bays_sub.sh $om0 $om1 $core &
done
wait
date

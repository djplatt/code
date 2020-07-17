#!/bin/bash
#
#PBS -l walltime=12:00:00,nodes=1:ppn=16
#PBS -q veryshort
date
hostname
q=$1
mkdir $HOME/zeros/${q}
for((t=0;t<220;t+=5))
do
   t1=`expr $t + 5`
   $HOME/code/l-func-hi/l-func $q $HOME/data/hurwitz_${t}_${t1}.dat $HOME/zeros/${q}/l-func_${t}_${t1}.dat
done


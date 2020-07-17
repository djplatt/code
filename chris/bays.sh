#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
zd="${HOME}/data/zeros"
n=0
for ((om=4377823041;om<=4377823059;om+=1))
do
   mkdir -p $HOME/data/chris/bays_437782/$om
   for t in 14 5000 26000 236000 446000
   do
       ${HOME}/code/chris/M_sum 200 ${zd}/zeros_${t}.dat ${om} 10000000 1.64e19 > ${HOME}/data/chris/bays_437782/${om}/${t}.out &
   done
   n=`expr $n + 1`
   if [ $n -ge 50 ]
   then
      wait
      n=0
   fi
done
wait
date

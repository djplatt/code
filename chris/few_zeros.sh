#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
zd="${HOME}/data/zeros"
n=0
for ((om=$om_start;om<$om_end;om+=1))
do
   mkdir $HOME/data/chris/few/$om
   for t in 14 5000 26000 236000 446000
   do
       ${HOME}/code/chris/M_sum 200 ${zd}/zeros_${t}.dat ${om} 1000 2.8363226783226e11 > ${HOME}/data/chris/few/${om}/${t}.out &
   done
   n=`expr $n + 1`
   if [ $n -ge 50 ]
   then
      wait
      n=0
   fi
done
wait
for ((om=$om_start;om<$om_end;om+=1))
do
   echo "$om $(grep "test on" ${HOME}/data/chris/few/$om/*.out | awk '//{print $19 " " $20 " " $21}' | $HOME/code/chris/arb_sum)"
done > $HOME/data/chris/few/few_${om_start}_${om_end}.out
date

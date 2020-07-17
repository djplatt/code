#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
ncores=16
for((om=$om_start;om<=$om_end;om+=2))
do
   mkdir -p $HOME/data/chris/billion/$om
   n=0
   for f in $(cat ~/data/zeros/first_billion.lst)
   do
       ${HOME}/code/chris/M_sum 200 $HOME/data/zeros/${f} ${om} 10000000000 3.599023e19 > ${HOME}/data/chris/billion/${om}/${n}.out &
       n=`expr $n + 1`
   done
   wait
done
for((om=$om_start;om<=$om_end;om++))
do
   echo "$om $(cat ${HOME}/data/chris/billion/$om/*.out | awk '//{print $19 " " $20 " " $21}' | ${HOME}/code/chris/arb_sum)" >> $HOME/data/chris/billion_${om_start}_${om_end}.lst
done
date

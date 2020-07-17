#!/bin/bash
#
#PBS -l walltime=10:00:00,nodes=1:ppn=16
#PBS -q veryshort
#
date
cd ${HOME}/code/yannick
rm -f ${list}.out
for f in $(cat $list)
do
   for g in $(cat $f)
   do
      $HOME/code/yannick/skewes_arb 200 ${HOME}/data/zeros/${g} 727951332982 1000000000 1.0e18 >> ${list}.out &
   done
   wait
done
date

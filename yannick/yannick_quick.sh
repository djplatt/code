#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
#
date
cd ${HOME}/code/yannick
for((om=727951332980;om<=727951332990;om++))
do
rm -f results/${om}_34046000.out
for g in 14 5000 26000 236000
 do
    $HOME/code/yannick/skewes_arb 200 ${HOME}/data/zeros/zeros_${g}.dat $om 1000000000 1.0e18 >> results/${om}_34046000.out &
 done
for ((g=446000;g<34046000;g+=2100000))
 do
    $HOME/code/yannick/skewes_arb 200 ${HOME}/data/zeros/zeros_${g}.dat $om 1000000000 1.0e18 >> results/${om}_34046000.out &
 done
wait
done
date

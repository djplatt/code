#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
cd $HOME/data/chris
for f in $(ls x*)
do
   $HOME/shells/skewes_chris1.sh $f &
done
wait
date

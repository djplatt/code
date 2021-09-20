#!/bin/bash
#
#PBS -l walltime=96:00:00,nodes=1:ppn=16
#PBS -q short
#
date
cd ${HOME}/code/harald/apr_20
for t in $(cat $list)
do
      $HOME/code/harald/apr_20/residues 150 ${HOME}/data/zeros/zeros_${t}.dat data/residues_${t}.out &
done
wait
date

#!/bin/bash
#
#PBS -l walltime=6:00:00,nodes=1:ppn=16
#PBS -q veryshort
#
#---------------------------------------------
# you would edit this section
# 
#
#---------------------------------------------
#   everything from here on is standard
#
date
cd /newhome/madjp/code/Q/runs
for f in $(ls x??)
do
   ./do_it.sh $f > ${f}.out &
done
wait

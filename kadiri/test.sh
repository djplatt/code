#!/bin/bash
#
#PBS -l walltime=0:10:00,nodes=1:ppn=8
#PBS -q testq
#
#---------------------------------------------
# you would edit this section
# 
#
#---------------------------------------------
#   everything from here on is standard
#
date
hostname
cd $HOME/code/kadiri
rm -f zeros1_sum_*.txt
n=0
for fn in $(ls $HOME/zeta_zeros/zeros_140*.dat)
do
   ./sum_rho $fn >> zeros1_sum_${n}.txt &
   n=`expr $n + 1`
   if [ $n -eq 8 ]
   then
      wait
      n=0
   fi
done
wait
grep igma zeros1_sum*.txt | awk '{print $3}' > zeros1_sum.txt
echo "[0.0,0.0]" >> zeros1_sum.txt
./sum_mpfi zeros1_sum.txt 300 > zeros1_sum.res
date

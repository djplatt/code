#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
#
date
cd $HOME/code/thomas/kacz
if [ -z $st ]
   then
   echo "st not specified"
   exit
fi
step=2500
step1=$(($step * 10))
st1=$(($st * 10))
for((n=0;n<16;n++))
do
   echo "t0=$st1;t1=$(($st1 + $step1));" > kacz_${n}.gp
   cat kacz.gp >> kacz_${n}.gp
   st1=$(($st1 + $step1))
done
for((n=0;n<16;n++))
do
   gp -q kacz_${n}.gp > kacz_${st}_${n}.out &
done
wait
for((n=0;n<16;n++))
do
   rm kacz_${n}.gp
done
date

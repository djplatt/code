#!/bin/bash
#PBS -P xh3
#PBS -l walltime=1:00:00,mem=32GB,ncpus=16
#PBS -q express
date
mkdir -p data/test
cd data/test
step=`expr $nits \* $wid`
for((n=0;n<$threads;n++))
do
   $HOME/code/arb_zeta/arb_zeta $prec $st $nits $wid > arb_zeta_${st}.out &
   st=`expr $st + $step`
done
wait
date

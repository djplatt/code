#!/bin/bash
#PBS -P qk9
#PBS -l walltime=24:00:00,mem=128GB,ncpus=16
#PBS -q normal
date
hostname
mkdir -p data/arb_zeta
cd data/arb_zeta
step=`expr $nits \* $wid`
for((n=0;n<16;n++))
do
   $HOME/code/arb_zeta/arb_zeta $prec $st $nits $wid > arb_zeta_${st}.out &
   st=`expr $st + $step`
done
wait
date

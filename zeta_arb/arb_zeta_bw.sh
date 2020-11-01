#!/bin/bash
#PBS -P qk9
#PBS -l walltime=18:00:00,mem=256GB,ncpus=28
#PBS -q normalbw
date
hostname
mkdir -p data/arb_zeta_bw
cd data/arb_zeta_bw
step=`expr $nits \* $wid`
for((n=0;n<28;n++))
do
   $HOME/code/arb_zeta/arb_zeta $prec $st $nits $wid > arb_zeta_${st}.out &
   st=`expr $st + $step`
done
wait
date

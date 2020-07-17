#!/bin/bash
#
# $1 = q
#
q=$1
num_files=32
prec=100
specfile=$HOME/temp/spec_$q.lst
echo $num_files > $specfile
for((n=0;n<$num_files;n++))
do
   ofile=$HOME/temp/f_even_${q}_${n}.trm
   ofile1=${ofile}.dat
   echo $ofile1 >> $specfile
   $HOME/code/ramare/f_hat_even_terms $prec $q $num_files $n $ofile
done
for f in $(ls $HOME/temp/f_even_${q}*.trm)
do
   $HOME/code/ramare/f_hat_even1.1 $f $f.dat
done
$HOME/code/ramare/f_even $specfile $HOME/temp/zeros_${q} 1 0
for f in $(ls $HOME/temp/zeros_${q}*)
do
   $HOME/code/ramare/low_upsam8-1.1.10 $f $f.8
done


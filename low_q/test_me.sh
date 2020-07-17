#!/bin/bash
q=$1
nfiles=32
resdir=$HOME/test/test_$q
mkdir $resdir
echo $nfiles > $resdir/even_spec.lst
echo $nfiles > $resdir/odd_spec.lst
for((n=0;n<$nfiles;n++))
do
   echo $resdir/f_even_${q}_${n}.trm.dat >> $resdir/even_spec.lst
   $HOME/code/low_q/f_hat_even_terms1.3-static 100 $q $nfiles $n $resdir/f_even_${q}_${n}.trm
   $HOME/code/low_q/f_hat_even1.6-static $resdir/f_even_${q}_${n}.trm $resdir/f_even_${q}_${n}.trm.dat
   rm $resdir/f_even_${q}_${n}.trm
   echo $resdir/f_odd_${q}_${n}.trm.dat >> $resdir/odd_spec.lst
   $HOME/code/low_q/f_hat_odd_terms1.3-static 100 $q $nfiles $n $resdir/f_odd_${q}_${n}.trm
   $HOME/code/low_q/f_hat_odd1.6-static $resdir/f_odd_${q}_${n}.trm $resdir/f_odd_${q}_${n}.trm.dat
   rm $resdir/f_odd_${q}_${n}.trm
done
$HOME/code/low_q/f_even1.3 $resdir/even_spec.lst zeros 1 0
rm $resdir/f_even_*_*.trm.dat
$HOME/code/low_q/f_odd1.3 $resdir/odd_spec.lst zeros 1 0
rm $resdir/f_odd_*_*.trm.dat
   for f in $(ls $resdir/zeros*)
     do
     $HOME/code/low_q/low_upsam8-1.1.10-static $f $f.8 > $f.8.log
     done
   for f in $(ls $resdir/zeros*.8)
     do
     if [ -s $f ]
         then
	 $HOME/code/low_q/low_upsammoreC32-static $f $f.32 > $f.32.log
     fi
     done
   for f in $(ls $resdir/zeros*.8.32)
     do
     if [ -s $f ]
         then
	 $HOME/code/low_q/low_upsammoreC128-static $f $f.128 > $f.128.log
     fi
     done
   for f in $(ls $resdir/zeros*.8.32.128)
     do
     if [ -s $f ]
         then
	 $HOME/code/low_q/low_upsammoreD512-static $f $f.512 > $f.512.log
     fi
     done
   for f in $(ls $resdir/zeros*.8.32.128.512)
     do
     if [ -s $f ]
         then
	 $HOME/code/zeros/upsamdouble-static $f $f.D > $f.D.log
     fi
     done


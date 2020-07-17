#!/bin/bash
cd $HOME/temp
qs=1000
nfiles="32"
for q in $qs
do
   sfn="spec_${q}.txt"
   echo $nfiles > $sfn
#N_last=`expr $N_BY_2 - $STEP`
   nc=0
   echo "Starting..."
   date
   echo "Running f_hat_terms..."
#for((n0=0;n0<$STEP;n0+=$STEP1))
#do
#   n1=`expr $n0 + $STEP1`
#   ofn="f_hat_${q}_${n0}_${n1}"
#   $HOME/code/low_q/f_hat_odd_terms 100 $q $eta $n0 $n1 $N ${ofn}.trm &
#done
#wait
   for((n0=0;n0<$nfiles;n0++))
     do
     ofn="f_hat_${q}_${n0}"
     echo ${ofn}.trm.dat >> $sfn
#     echo "$HOME/code/low_q/f_hat_even_terms1.2 100 $q $nfiles $n0 ${ofn}.trm"
     $HOME/code/low_q/temp 100 $q $nfiles $n0 ${ofn}.trm
   done
   date
   echo "Running f_hat..."
   nc=0
   for ofn in $(ls *.trm)
      do
      $HOME/code/low_q/f_hat_even1.1 ${ofn} $HOME/data/facs_file.dat ${ofn}.dat
   done
#   rm *.trm
   date
   echo "Running f..."
   $HOME/code/low_q/f_even1.2 $sfn zeros0 1 0
#   rm *.trm.dat
   date
#   cp zeros*.dat $RESDIR/.
#   date
#   rm *
#   cd ..
#   rmdir $WORKDIR
#   exit
#
#
   echo "Running low_upsam..."
   for fn in $(ls zeros*${q}*.dat)
   do
      $HOME/code/low_q/low_upsam8-1.1.10 ${fn} ${fn}.8 > ${fn}.8.log
   done
#   ls -al
   date
   echo "Running low_upsammoreC32..."
   for fn in $(ls -al zeros*${q}*.dat.8 | grep " 0 [JFMASOND]" -v | awk '// {print $9}')
   do
      $HOME/code/low_q/low_upsammoreC32 ${fn} ${fn}.32 > ${fn}.32.log
   done
   date
   echo "Running low_upsammoreC128..."
   for fn in $(ls -al zeros*${q}*.dat.8.32 | grep " 0 [JFMASOND]" -v | awk '// {print $9}')
   do
      $HOME/code/low_q/low_upsammoreC128 ${fn} ${fn}.128 > ${fn}.128.log
   done
   date
   echo "Running low_upsammoreC512..."
   for fn in $(ls -al zeros*${q}*.dat.8.32.128 | grep " 0 [JFMASOND]" -v | awk '// {print $9}')
   do
      $HOME/code/low_q/low_upsammoreD512 ${fn} ${fn}.512 > ${fn}.512.log
   done
   date
   echo "Running upsamdouble..."
   for fn in $(ls -al zeros*${q}*.dat.8.32.128.512 | grep " 0 [JFMASOND]" -v | awk '// {print $9}')
      do
      $HOME/code/zeros/upsamdouble ${fn} facs_file.dat ${fn}.D > ${fn}.D.log
   done
   echo "Finished ${q}."
   date
done

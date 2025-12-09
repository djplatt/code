#!/bin/bash
fstub="res_${1}_${2}"
awk -v st1=$1 -v st2=$2 '//{if(($2>=st1)&&($2<st2)){printf "%d %d %d\n",st,st+10,int($5+0.5);st=st+10}} BEGIN {st=st1}' results.sum > ${fstub}.exp
#
awk '/ero/{if($16<1e-40){print substr($9,2,10) " + " substr($11,1,12) "i"}}' ${fstub}.txt | sort -n -k 2 | sort -u | sort -k 3 -n > ${fstub}.uni
#
rm -f $fstub.act
for((st=$1;st<$2;st+=10))
do
    echo "$st $(($st + 10)) $(awk -v st=$st '//{if((int($3)>=st)&&(int($3)<st+10)){count=count+1}} BEFIN {count=0} END {print count}' ${fstub}.uni)" >> $fstub.act
done

#!/bin/bash
for((n=1;n<31;n++))
do
    nn=$((2**${n}))
    echo "$n " $(awk -v N=${nn} '//{if($1<=N){count+=1}} BEGIN {count=0} END {print count}' 2_30.out)
done

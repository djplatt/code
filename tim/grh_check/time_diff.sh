#!/bin/bash
a=$(head -n 1 $1 | awk '//{print $4}')
b=$(tail -n 1 $1 | awk '//{print $4}')
as=$(date -d $a +%s)
bs=$(date -d $b +%s)
let delta=$bs-$as
if [ $delta -lt 0 ]
then
   let delta=$delta+84600
fi
echo "$delta"

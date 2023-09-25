#!/bin/bash
date
mkdir -p /user/work/madjp/watkins
cd /user/home/madjp/code/tim/watkins
for((start=5;start<300000000;start+=10714286))
do
  end=$(($start+10714285))
  ./watkins_seg_int_double $start $end > /user/work/madjp/watkins/${start}_${end}.out &
done
wait
date

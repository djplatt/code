#!/bin/bash
date
for((c=0;c<8;c++))
do
    batch_me1.sh $c &
done
wait
date

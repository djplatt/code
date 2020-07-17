#!/bin/sh
counter=0
while [ $counter -lt 228 ]
do
   ./test_local 100000000 &
   counter=`expr $counter + 1`
done
wait


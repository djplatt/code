#!/bin/bash
#
for f in $(cat $1)
do
   $HOME/code/chris/skewes 200 $f 727951335420 1000000000 1.0e16 >> $1.out
done

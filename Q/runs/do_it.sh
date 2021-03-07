#!/bin/bash
for f in $(cat $1)
do
   ../Q /newhome/madjp/data/zeros/${f} | tail -n 3
done

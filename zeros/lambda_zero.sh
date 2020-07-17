#!/bin/bash
#
#
for qi in "\"83183 8593\" \"86297 35148\""
do
  echo $qi
  q=${qi:0:5}
  ~/work/current-21-11/code/zeros/lambda_zero 100 ~/work/current-21-11/data/facs_file.dat ${qi} > ${q}.log
done

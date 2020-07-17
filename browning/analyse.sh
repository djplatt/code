#!/bin/bash
for((B=1000;B<=256000;B+=B))
do
   echo "$B: $(awk -v B=$B '//{if($2*$2<=B*B && $3*$3<=B*B && $4*$4<=B*B && $5*$5<=B*B) {print $2 " " $3 " " $4 " " $5}}' $1 | wc -l)"
done
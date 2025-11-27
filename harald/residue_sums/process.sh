#!/bin/bash
#
# strip out the 6 sums in arb_load_file format
# output by residue_sums_nofft
# assumes the output from each residue file has been cat'd together
grep --no-group-separator -A 6 "We processed" $1 | grep "We processed" -v > temp
split -n r/6 temp $1
rm -rf temp
for f in ${1}??
do
    grep "0 0 0 0" -v $f > $f.out
    echo "0 0 0 0" >> $f.out
    cat $f.out | ./arb_load_sum 150
done



#!/bin/bash
a=$(head -n 1 $1)
b=$(tail -n 1 $1)
as=$(date -d "$a" +%s)
bs=$(date -d "$b" +%s)
#echo "$a $b $as $bs"
delta="$(($bs-$as))"
echo $delta


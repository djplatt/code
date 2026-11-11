#!/bin/bash
c=$1
st=$((1000000 * ($c - 1)))
en=$(($st + 999999))
for((n=0;n<133;n++))
do
    echo "../sieve_only $st $en > $st_$en.out"
    st=$(($st + 8000000))
    en=$(($en + 8000000))
done

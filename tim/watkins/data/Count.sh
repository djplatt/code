#!/bin/bash
grep -h "There are" */*.out | sort -n -k 8 | awk '//{if(cn!=$8) {print cn " " count;cn=$8;count=$3} else {count+=$3}} BEGIN {count=0;cn="foobar"} END {print cn " " count}'

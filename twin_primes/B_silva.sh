#!/bin/bash
for((n=9;n<=15;n++)); do ./B_silva 2d${n}.dat 200; done | awk '//{sum1+=$2;sum+=$5} END {printf "%10.8e %10.8e\n",sum1+1.8065924,sum+1.8065925}'

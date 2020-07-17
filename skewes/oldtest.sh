#!/bin/bash
#
#PBS -l walltime=12:00:00,nodes=1:ppn=8
#PBS -q veryshort
#
#---------------------------------------------
# you would edit this section
# 
#
#---------------------------------------------
#   everything from here on is standard
#
#  get name for scratch directory and create it
#
echo "job started on"
hostname
date
$HOME/code/skewes/check_li1.2 10000000000000 10100000000000 346065536839 349405678546 $cs &
$HOME/code/skewes/check_li1.2 10000000000000 10100000000000 346065536839 349405678546 $cs &
$HOME/code/skewes/check_li1.2 10000000000000 10100000000000 346065536839 349405678546 $cs &
$HOME/code/skewes/check_li1.2 10000000000000 10100000000000 346065536839 349405678546 $cs &
$HOME/code/skewes/check_li1.2 10000000000000 10100000000000 346065536839 349405678546 $cs &
$HOME/code/skewes/check_li1.2 10000000000000 10100000000000 346065536839 349405678546 $cs &
$HOME/code/skewes/check_li1.2 10000000000000 10100000000000 346065536839 349405678546 $cs &
$HOME/code/skewes/check_li1.2 10000000000000 10100000000000 346065536839 349405678546 $cs &
wait
date

#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
time $HOME/primesieve-5.4.1/primesieve -c2 -t32 100000000000000000 100010000000000000
date

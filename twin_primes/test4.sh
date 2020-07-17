#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
cd $HOME/code/twin_primes
rm -f results4.dat
time ./twin_primes 1e17 10000128e10 1e10 >> results4.dat
date

#!/bin/bash
#
#PBS -l walltime=1:00:00,nodes=1:ppn=16
#PBS -q testq
date
hostname
time $HOME/code/farey/farey_threads_int 1000000 16
date

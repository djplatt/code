#!/bin/bash
#
#PBS -l walltime=25:00:00,nodes=1:ppn=16
#PBS -q short
date
hostname
time $HOME/code/farey/farey_threads_int $N 16
date

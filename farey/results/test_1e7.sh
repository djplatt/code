#!/bin/bash
#
#PBS -l walltime=10:00:00,nodes=1:ppn=16
#PBS -q veryshort
date
hostname
time $HOME/code/farey/farey_threads 10000000 16
date

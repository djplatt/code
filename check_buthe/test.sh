#!/bin/bash
#
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=1
#PBS -q short
#
date
hostname
cd ${HOME}/code/check_buthe
./zeta-zeros btmp djp_zeros > check.log
date

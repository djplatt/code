#!/bin/bash
for f in x??
do
   qsub -joe -v dfile=$f do_t.sh
done
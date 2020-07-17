#!/bin/bash
#
#PBS -l walltime=1:00:00
#
#---------------------------------------------
# you would edit this section
#
setenv MYDIR "${HOME}/data"
setenv MYEX "${HOME}/code/l-func-mpfi 300 50 30 30 200 4096 z_file_five.dat foo.dat 15 1"
#
#---------------------------------------------
#   everything from here on is standard
#
#  get name for scratch directory and create it
#
setenv JOBNO "`echo $PBS_JOBID | sed s/.bluequeue1.cvos.cluster//`"
setenv WORKDIR "/local/${PBS_O_LOGNAME}.${JOBNO}"
mkdir $WORKDIR
#
setenv MYDIR "${HOME}/KRUNS/GNU/run02"
setenv MYEX "${HOME}/KSPACE/GNU/bin/KSPACE"
#

#  go to scratch directory and copy files to it
#
cd $WORKDIR
cp ${MYDIR}/z_file_five.dat .
#
#  store name of node and scratch directory in file
#  rundat in file directory
#
hostname > ${MYDIR}/rundat
pwd >> ${MYDIR}/rundat
#
#   run program
#
$MYEX
#
#  copy files back and delete scratch directory
#
cd ${MYDIR}
mv -f ${WORKDIR}/* .
rm -fr $WORKDIR

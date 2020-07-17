#!/bin/bash
CRDIR="$HOME/crlibm-1.0beta4/"
#echo $CRDIR
#echo $MPFI_LIB
#echo $MPFI_INCLUDE
PARI_INCLUDE=/usr/local/apps/pari-2.4.3.alpha/include/pari
PARI_LIB=/usr/local/apps/pari-2.4.3.alpha/lib
g++ -o$HOME/code/skewes/check_li1.2 $HOME/code/skewes/check_li1.2.cpp $HOME/code/primesieve/out/*.o -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -I$MPFI_INCLUDE -L$MPFI_LIB -I$FFTW3_INCLUDE -I$PARI_INCLUDE -L$PARI_LIB -L$FFTW3_LIB -lmpfi -lcrlibm -lrt -lfftw3 -lpari

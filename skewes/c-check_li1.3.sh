#!/bin/bash
CRDIR="$HOME/crlibm-1.0beta4/"
#echo $CRDIR
#echo $MPFI_LIB
#echo $MPFI_INCLUDE
#PARI_INCLUDE=/usr/local/apps/pari-2.4.3.alpha/include/pari
#PARI_LIB=/usr/local/apps/pari-2.4.3.alpha/lib
g++ -o$HOME/code/skewes/check_li1.3 $HOME/code/skewes/check_li1.3.cpp $HOME/code/primesieve/out/*.o -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -I$MPFI_INCLUDE -L$MPFI_LIB -lmpfi -lcrlibm

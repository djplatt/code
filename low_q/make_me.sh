#!/bin/bash
#
# Poor man's makefile
#
# Don't go beyond -O1, things tend to break!
# C files need MPFI (and therefore MPFR and GMP)
gcc f_hat_even_terms1.3.c -of_hat_even_terms1.3-static -static -O1 -finline-functions -fomit-frame-pointer -I$MPFI_INCLUDE -L$MPFI_LIB -lmpfi -lmpfr -lgmp -lm
gcc f_hat_odd_terms1.3.c -of_hat_odd_terms1.3-static -static -O1 -finline-functions -fomit-frame-pointer -I$MPFI_INCLUDE -L$MPFI_LIB -lmpfi -lmpfr -lgmp -lm
gcc f_hat_even_terms1.3.c -of_hat_even_terms1.3 -O1 -finline-functions -fomit-frame-pointer -I$MPFI_INCLUDE -L$MPFI_LIB -lmpfi -lm
gcc f_hat_odd_terms1.3.c -of_hat_odd_terms1.3 -O1 -finline-functions -fomit-frame-pointer -I$MPFI_INCLUDE -L$MPFI_LIB -lmpfi -lm

# C++ files need crlibm
# this is where I keep crlibm
CRDIR="$HOME/crlibm-1.0beta4/"
g++ f_hat_even1.6.cpp -of_hat_even1.6 -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ f_even1.3.cpp -of_even1.3 -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ f_hat_odd1.6.cpp -of_hat_odd1.6 -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ f_odd1.3.cpp -of_odd1.3 -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ low_upsam8-1.1.10.cpp -olow_upsam8-1.1.10 -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ low_upsammoreC32.cpp -olow_upsammoreC32 -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ low_upsammoreC128.cpp -olow_upsammoreC128 -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ low_upsammoreD512.cpp -olow_upsammoreD512 -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ ../zeros/upsamdouble.cpp -o../zeros/upsamdouble -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ f1.3.cpp -of1.3 -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ f_hat_split1.2.cpp -of_hat_split1.2 -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
#
g++ f_hat_even1.6.cpp -of_hat_even1.6-static -static -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ f_even1.3.cpp -of_even1.3-static -static -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ f_hat_odd1.6.cpp -of_hat_odd1.6-static -static -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ f_odd1.3.cpp -of_odd1.3-static -static -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ low_upsam8-1.1.10.cpp -olow_upsam8-1.1.10-static -static -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ low_upsammoreC32.cpp -olow_upsammoreC32-static -static -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ low_upsammoreC128.cpp -olow_upsammoreC128-static -static -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ low_upsammoreD512.cpp -olow_upsammoreD512-static -static -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ ../zeros/upsamdouble.cpp -o../zeros/upsamdouble-static -static -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ f1.3.cpp -of1.3-static -static -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm
g++ f_hat_split1.2.cpp -of_hat_split1.2-static -static -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lm


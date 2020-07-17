gcc -O2 -o sum_rho sum_rho.c -I ~/mpfi-1.5.1/src -I ~/mpfr-3.1.2/src -L ~/mpfi-1.5.1/src/.libs -lmpfi -lm
gcc -O2 -o sum_mpfi sum_mpfi.c -I ~/mpfi-1.5.1/src -I ~/mpfr-3.1.2/src -L ~/mpfi-1.5.1/src/.libs -lmpfi -lm

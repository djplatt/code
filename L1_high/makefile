all:	L1_high print_L1_high

#L1_high: L1_high.cpp ../includes/mpfi_c_hurwitz.h
#	g++ -O2 -o L1_high L1_high.cpp ../characters-master/mpfi-fft.o -I ${FFTWINCLUDE} -I /panfs/panasas01/math/madjp/arb-master/ -I /panfs/panasas01/math/madjp/flint2-trunk -I /panfs/panasas01/math/madjp/mpfi-1.5.1/src -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/gmp-5.1.3 -L /panfs/panasas01/math/madjp/arb-master/ -larb -L /panfs/panasas01/math/madjp/flint2-trunk -lflint -L /panfs/panasas01/math/madjp/mpfi-1.5.1/src/.libs/ -lmpfi -L /panfs/panasas01/math/madjp/mpfr-3.1.2/src/.libs/ -lmpfr -L /panfs/panasas01/math/madjp/gmp-5.1.3/.libs/ -lgmp -L ${FFTWDIR} -lfftw3

print_L1_high: print_L1_high.cpp
	g++ -O2 -o print_L1_high print_L1_high.cpp -lgmp

L1_high:	L1_high.cpp
	g++ -O2 L1_high.cpp ../characters-master/mpfi-fft.o -o L1_high -I ${FFTWINCLUDE} -I /panfs/panasas01/math/madjp/mpfi-1.5.1/src -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.8.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.8.1 -larb -L /panfs/panasas01/math/madjp/mpfi-1.5.1/src/.libs/ -lmpfi -L ${FFTWDIR} -lfftw3

#print_L1_high:	print_L1_high.cpp
#	g++ -O2 print_L1_high.cpp ../characters-master/mpfi-fft.o -o print_L1_high -I ${FFTWINCLUDE} -I /panfs/panasas01/math/madjp/mpfi-1.5.1/src -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.8.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.8.1 -larb -L /panfs/panasas01/math/madjp/mpfi-1.5.1/src/.libs/ -lmpfi -L ${FFTWDIR} -lfftw3
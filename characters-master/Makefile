all: chimax dft-test test-values test-primitive

chimax: chimax.cc characters_djp.h slint.h mpfi-fft.o
	g++ -O3 -o chimax chimax.cc -pthread mpfi-fft.o -I${MPFI_INCLUDE} -lfftw3 -lmpfi -lmpfr
dft-test: dft-test.cc characters_djp.h slint.h mpfi-fft.o
	g++ -O3 -o dft-test dft-test.cc mpfi-fft.o -lfftw3 -lmpfi -lmpfr
test-primitive: characters_djp.h slint.h test-primitive.cc mpfi-fft.o
	g++ -O3 -o test-primitive test-primitive.cc mpfi-fft.o -lfftw3 -lmpfi -lmpfr
test-values: characters_djp.h slint.h test-values.cc mpfi-fft.o
	g++ -O2 -o test-values test-values.cc mpfi-fft.o -lfftw3 -lmpfi -lmpfr
test-mpfi-fft: mpfi-fft.o test-mpfi-fft.o
	gcc -o test-mpfi-fft mpfi-fft.o test-mpfi-fft.o -lmpfi -lmpfr
mpfi-fft.o: mpfi_fft.h mpfi_fft.c
	gcc -c -o mpfi-fft.o mpfi_fft.c -O2 -I /panfs/panasas01/math/madjp/mpfi-1.5.1/src 
test-mpfi-fft.o: test-mpfi-fft.c mpfi_fft.h
	gcc -c -o test-mpfi-fft.o test-mpfi-fft.c -O2

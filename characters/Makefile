all: test-mpfi-fft mpfi-fft.o

ROOT=.
INCLUDE=${ROOT}
MPFI_INCLUDE=-I ${HOME}/mpfi-1.5.1/src -I ${HOME}/mpfr-3.1.2/src
MPFI_LIB=-L ${HOME}/mpfi-1.5.1/src/.libs

test-mpfi-fft: mpfi-fft.o test-mpfi-fft.o
	gcc -o test-mpfi-fft mpfi-fft.o test-mpfi-fft.o ${MPFI_LIB} -lmpfi
mpfi-fft.o: $(INCLUDE)/mpfi_fft.h mpfi_fft.c
	gcc -c -o mpfi-fft.o mpfi_fft.c -O2 ${MPFI_INCLUDE}
test-mpfi-fft.o: test-mpfi-fft.c $(INCLUDE)/mpfi_fft.h
	gcc -c -o test-mpfi-fft.o test-mpfi-fft.c -O2 ${MPFI_INCLUDE}

clean:
	-rm test-mpfi-fft
	-rm mpfi-fft.o
	-rm test-mpfi-fft.o

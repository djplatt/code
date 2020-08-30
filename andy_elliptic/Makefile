LIB=/data/local/lib
INC=/data/local/include

all: plot2 g2_coeffs elliptic_coeffs

elliptic_coeffs: elliptic_coeffs.c bessel_arb.o fft_acb.o pari.c
	gcc -O -I$(INC) -L$(LIB) -o elliptic_coeffs elliptic_coeffs.c bessel_arb.o fft_acb.o -lgmp -lflint -larb -lm -lsmalljac -lff_poly -lpari

g2_coeffs: g2_coeffs.c bessel_arb.o fft_acb.o pari.c
	gcc -O -I$(INC) -L$(LIB) -o g2_coeffs g2_coeffs.c bessel_arb.o fft_acb.o -lgmp -lflint -larb -lm -lsmalljac -lff_poly -lpari


plot2: plot2.c bessel_arb.o fft_acb.o pari.c
	gcc -O -I$(INC) -L$(LIB) -o plot2 plot2.c bessel_arb.o fft_acb.o -lgmp -lflint -larb -lm -lsmalljac -lff_poly -lpari

bessel_arb.o: bessel_arb.c wdata.c
	gcc -O -I$(INC) -c bessel_arb.c

fft_acb.o: fft_acb.c
	gcc -O -I$(INC) -c fft_acb.c

clean:
	rm -f *.o plot2 g2_coeffs

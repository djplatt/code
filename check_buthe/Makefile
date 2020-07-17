CFLAGS=-O2 -finline-functions -fomit-frame-pointer -I ~/gmp-5.1.3 -I ~/mpfr-3.1.2/src -I ~/arb-2.8.1/ -I ~/flint2-trunk
LDFLAGS=-L ~/mpfr-3.1.2/src/.libs -L ~/arb-2.8.1/ -larb
CC = gcc -g $(CFLAGS) $(LDFLAGS)
CPP = g++ $(CFLAGS) $(LDFLAGS)
all: zeta-zeros check_one list_buthe_zeros zeta-zeros-265

zeta-zeros: zero_stream.c zeta-zeros.c if.o djp_zeros.c
	$(CC) zeta-zeros.c zero_stream.c if.o -larb -o zeta-zeros

zeta-zeros-265: zero_stream.c zeta-zeros-265.c if.o djp_zeros.c
	$(CC) zeta-zeros-265.c zero_stream.c if.o -larb -o zeta-zeros-265

list_buthe_zeros: zero_stream.c list_buthe_zeros.c if.o djp_zeros.c
	$(CC) list_buthe_zeros.c zero_stream.c if.o -larb -o list_buthe_zeros

check_one: check_one.cpp
	$(CPP) check_one.cpp -o check_one

clear:
	rm zeta-zeros if.o

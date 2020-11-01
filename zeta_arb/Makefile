all: arb_zeta_local

turing.o: turing.c turing.h parameters.h
	gcc -c -Wall turing.c -O2 -finline-functions -fomit-frame-pointer

inter.o: inter.c inter.h parameters.h
	gcc -c -Wall inter.c -O2 -finline-functions -fomit-frame-pointer

win_zeta.o: win_zeta.c win_zeta.h parameters.h
	gcc -c -Wall win_zeta.c -O2 -finline-functions -fomit-frame-pointer

arb_zeta.o: arb_zeta.c parameters.h
	gcc -c -Wall arb_zeta.c -O2 -finline-functions -fomit-frame-pointer

arb_zeta_local: arb_zeta.o win_zeta.o inter.o turing.o
	gcc -o arb_zeta arb_zeta.o turing.o inter.o win_zeta.o -l arb

clean:
	rm -f arb_zeta *.o
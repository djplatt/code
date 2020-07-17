INC=-I ${HOME}/code/includes -I ${HOME}/primesieve-5.4.1/include -I ${HOME}/arb-2.17.0 -I ${HOME}/flint-2.5.2/ -I ${HOME}/mpfr-3.1.2/src -I ${HOME}/gmp-5.1.3
LIB=-L ${HOME}/primesieve-5.4.1/.libs -L ${HOME}/arb-2.17.0

all:	zeta_prime

zeta_prime:	zeta_prime.cpp
	g++ -O2 zeta_prime.cpp -o zeta_prime ${INC} ${LIB} -larb


clean:
	rm -f zeta_prime

INC=-I ${HOME}/primesieve-5.4.1/include -I ${HOME}/arb-2.15.1 -I ${HOME}/flint-2.5.2/ -I ${HOME}/mpfr-3.1.2/src -I ${HOME}/gmp-5.1.3
LIB=-L ${HOME}/primesieve-5.4.1/.libs -L ${HOME}/arb-2.15.1

all:	saouter

saouter:	saouter.cpp
	g++ -O2 saouter.cpp -o saouter ${INC} ${LIB} -larb -lprimesieve

clean:
	rm -f saouter

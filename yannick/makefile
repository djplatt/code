INC=-I ${HOME}/code/includes -I ${HOME}/primesieve-5.4.1/include -I ${HOME}/arb-2.15.1 -I ${HOME}/flint-2.5.2/ -I ${HOME}/mpfr-3.1.2/src -I ${HOME}/gmp-5.1.3
LIB=-L ${HOME}/primesieve-5.4.1/.libs -L ${HOME}/arb-2.15.1

all:	skewes_arb arb_sum

skewes_arb:	skewes_arb.cpp
	g++ -O2 skewes_arb.cpp -o skewes_arb ${INC} ${LIB} -larb -lprimesieve

arb_sum:	arb_sum.cpp
	g++ -O2 arb_sum.cpp -o arb_sum ${INC} ${LIB} -larb -lprimesieve

clean:
	rm -f skewes_arb arb_sum

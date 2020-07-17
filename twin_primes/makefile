INC=${HOME}/primesieve-5.4.1/include
LIB=${HOME}/primesieve-5.4.1/.libs

all:	twin_primes twin_primes_serial H B B_silva B_int

twin_primes:	twin_primes.cpp
	g++ -O2 twin_primes.cpp -o twin_primes -I ${INC} -L ${LIB} -lprimesieve

twin_primes_serial:	twin_primes_serial.cpp
	g++ -O2 twin_primes_serial.cpp -o twin_primes_serial -I ${INC} -L ${LIB}  -lprimesieve

H:	H.cpp
	g++ -O2 H.cpp -o H -I ${INC} -L ${LIB} -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.8.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.8.1 -larb -lprimesieve

B:	B.cpp
	g++ -O2 B.cpp -o B -I ${INC} -L ${LIB} -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.8.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.8.1 -larb -lprimesieve

B_silva:	B_silva.cpp
	g++ -O2 B_silva.cpp -o B_silva -I ${INC} -L ${LIB} -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.8.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.8.1 -larb

B_int:	B_int_fixed.cpp consts.h
	g++ -O2 B_int_fixed.cpp -o B_int -I ${INC} -L ${LIB} -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.8.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.8.1 -larb

clean:
	rm -f twin_primes
	rm -f twin_primes_serial
	rm -f H
	rm -f B
	rm -f B_silva
	rm -f B_int

PS_INC=-I /panfs/panasas01/math/madjp/primesieve-5.4.1/include
PS_LIB=-L /panfs/panasas01/math/madjp/primesieve-5.4.1/.libs
CR_INC=-I /panfs/panasas01/math/madjp/crlibm-1.0beta4
CR_LIB=-L /panfs/panasas01/math/madjp/crlibm-1.0beta4
ARB_INC=-I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.8.1 -I /panfs/panasas01/math/madjp/flint2-trunk
ARB_LIB=-L /panfs/panasas01/math/madjp/arb-2.8.1

all: int_double double long_double int_double_smart arb

arb: arb.cpp
	g++ arb.cpp -o arb -O1 -finline-functions -fomit-frame-pointer ${ARB_INC} ${PS_INC} ${ARB_LIB} -larb ${PS_LIB} -lprimesieve

int_double: int_double.cpp int_double.h
	g++ int_double.cpp -o int_double -O1 -frounding-math -finline-functions -fomit-frame-pointer ${PS_INC} ${CR_INC} ${PS_LIB} -lprimesieve ${CR_LIB} -lcrlibm

int_double_smart: int_double_smart.cpp int_double.h
	g++ int_double_smart.cpp -o int_double_smart -O1 -frounding-math -finline-functions -fomit-frame-pointer ${PS_INC} ${CR_INC} ${PS_LIB} -lprimesieve ${CR_LIB} -lcrlibm

double: double.cpp
	g++ double.cpp -o double -O1 -finline-functions -fomit-frame-pointer ${PS_INC} ${PS_LIB} -lprimesieve

long_double: long_double.cpp
	g++ long_double.cpp -o long_double -O1 -finline-functions -fomit-frame-pointer ${PS_INC} ${PS_LIB} -lprimesieve

clean:
	rm -f arb
	rm -f int_double
	rm -f double
	rm -f long_double
	rm -f int_double_smart

all: l_func l_func_old

l_func:	main8.0.cpp structs8.0.h error8.0.h acb_fft.h
	g++ main8.0.cpp -o l-func -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-master -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-master -larb 

l_func_old: main7.0.cpp structs7.0.h error7.0.h acb_fft.h
	g++ main7.0.cpp -o l-func-old -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-master -I /panfs/panasas01/math/madjp/flint2-trunk -I /panfs/panasas01/math/madjp/primesieve-5.4.1/include -I /panfs/panasas01/math/madjp/characters-master -L /panfs/panasas01/math/madjp/arb-master -larb -L /panfs/panasas01/math/madjp/primesieve-5.4.1/.libs -lprimesieve

test:	test.cpp acb_fft.h
	g++ test.cpp -o test -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-master -I /panfs/panasas01/math/madjp/flint2-trunk -I /panfs/panasas01/math/madjp/primesieve-5.4.1/include -I /panfs/panasas01/math/madjp/characters-master -L /panfs/panasas01/math/madjp/arb-master -larb -L /panfs/panasas01/math/madjp/primesieve-5.4.1/.libs -lprimesieve

all: chacha_riemann test_l_jet chacha2

chacha_riemann: chacha_riemann.cpp
	g++ -std=gnu++0x chacha_riemann.cpp -o chacha_riemann -O2 -finline-functions -fomit-frame-pointer -I ${HOME}/code/includes -I ${FFTWINCLUDE} -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1/ -I /panfs/panasas01/math/madjp/flint-2.5.2 -L/panfs/panasas01/math/madjp/gmp-5.1.3/.libs -L/panfs/panasas01/math/madjp/mpfr-3.1.2/src/.libs -L/panfs/panasas01/math/madjp/flint-2.5.2 -L /panfs/panasas01/math/madjp/arb-2.15.1/ -larb -lflint -lmpfr -lgmp -lm -lpthread

test_l_jet: test_l_jet.cpp
	g++ -std=gnu++0x test_l_jet.cpp -o test_l_jet -O2 -finline-functions -fomit-frame-pointer -I ${HOME}/code/includes -I ${FFTWINCLUDE} -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1/ -I /panfs/panasas01/math/madjp/flint-2.5.2 -L/panfs/panasas01/math/madjp/gmp-5.1.3/.libs -L/panfs/panasas01/math/madjp/mpfr-3.1.2/src/.libs -L/panfs/panasas01/math/madjp/flint-2.5.2 -L /panfs/panasas01/math/madjp/arb-2.15.1/ -larb -lflint -lmpfr -lgmp -lm

chacha2: chacha2.cpp
	g++ -std=gnu++0x chacha2.cpp -o chacha2 -O2 -finline-functions -fomit-frame-pointer -I ${HOME}/code/includes -I ${FFTWINCLUDE} -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1/ -I /panfs/panasas01/math/madjp/flint-2.5.2 -L/panfs/panasas01/math/madjp/gmp-5.1.3/.libs -L/panfs/panasas01/math/madjp/mpfr-3.1.2/src/.libs -L/panfs/panasas01/math/madjp/flint-2.5.2 -L /panfs/panasas01/math/madjp/arb-2.15.1/ -larb -lflint -lmpfr -lgmp -lm -lpthread

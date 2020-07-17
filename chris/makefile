all: skewes numint skewes_int skewes_arb arb_test mpfi_test M_sum arb_sum

arb_sum: arb_sum.cpp ../includes/arb_sum.h
	g++ arb_sum.cpp -o arb_sum -O1 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1 -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb 

skewes: skewes.c
	gcc -O2 -o skewes skewes.c -I ${HOME}/mpfi-1.5.1/src -L ${HOME}/mpfi-1.5.1/src/.libs -lmpfi -lm

mpfi_test: mpfi_test.cpp
	g++ -O2 -o mpfi_test mpfi_test.cpp -I ${HOME}/mpfi-1.5.1/src -L ${HOME}/mpfi-1.5.1/src/.libs -lmpfi -lm

numint: numint.cpp ../trace/quad.h
	g++ numint.cpp -o numint -O1 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1 -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb 

skewes_int: skewes_int.cpp ../trace/quad.h
	g++ skewes_int.cpp -o skewes_int -O1 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1 -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb 

skewes_arb: skewes_arb.cpp
	g++ skewes_arb.cpp -o skewes_arb -O1 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1 -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb 

arb_test: arb_test.cpp
	g++ arb_test.cpp -o arb_test -O1 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1 -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb 

M_sum: M_sum.cpp
	g++ M_sum.cpp -o M_sum -O1 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1 -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb 

rect_sum: rect_sum.cpp
	g++ rect_sum.cpp -o rect_sum -O1 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1 -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb 


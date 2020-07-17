all: arb_zeta_rj arb_zeta_bc

arb_zeta_bc: arb_zeta.cpp parameters.h
	g++ arb_zeta.cpp -o arb_zeta -O2 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.8.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.8.1 -larb

arb_zeta_rj: arb_zeta.cpp parameters.h
	echo "Not making for Raijin."

#	g++ arb_zeta.cpp -o arb_zeta -O2 -finline-functions -fomit-frame-pointer -I /apps/flint/2.5.2/include/flint/ -larb

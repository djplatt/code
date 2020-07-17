all: check_ram

check_ram: check_ram.cpp
	g++ check_ram.cpp -o check_ram -O2 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.8.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.8.1 -larb
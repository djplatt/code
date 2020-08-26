all: ex1 ex2


ex1: ex1.c
	gcc ex1.c -o ex1 -O2 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb

ex2: ex2.c
	gcc ex2.c -o ex2 -O2 -finline-functions -fomit-frame-pointer -I /panfs/panasas01/math/madjp/gmp-5.1.3 -I /panfs/panasas01/math/madjp/mpfr-3.1.2/src -I /panfs/panasas01/math/madjp/arb-2.15.1/ -I /panfs/panasas01/math/madjp/flint2-trunk -L /panfs/panasas01/math/madjp/arb-2.15.1 -larb

clean:
	rm -f ex2 ex1
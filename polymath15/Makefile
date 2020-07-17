all: H

H: H.c
	gcc H.c -o H -O1 -finline-functions -fomit-frame-pointer -I/panfs/panasas01/math/madjp/arb-2.13.0 -I/panfs/panasas01/math/madjp/gmp-5.1.3/ -I/panfs/panasas01/math/madjp/mpfr-3.1.2//src -I/panfs/panasas01/math/madjp/flint2-trunk -L/panfs/panasas01/math/madjp/arb-2.13.0 -L/panfs/panasas01/math/madjp/gmp-5.1.3//.libs -L/panfs/panasas01/math/madjp/mpfr-3.1.2//src/.libs -L/panfs/panasas01/math/madjp/flint2-trunk -larb -lflint -lmpfr -lgmp -lm -lpthread

clean:
	rm -f H

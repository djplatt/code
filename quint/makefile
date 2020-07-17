all: quintuples8 quintuples7 quintuples6 ineq

quintuples8:	quintuples8.c
	gcc quintuples8.c -o quintuples8 -O2 -I${HOME}/mpfi-1.5.1/src -I${HOME}/mpfr-3.1.2/src -L${HOME}/mpfi-1.5.1/src/.libs -lmpfi

quintuples7:	quintuples7.c
	gcc quintuples7.c -o quintuples7 -O2 -I${HOME}/mpfi-1.5.1/src -I${HOME}/mpfr-3.1.2/src -L${HOME}/mpfi-1.5.1/src/.libs -lmpfi


quintuples6:	quintuples6.c
	gcc quintuples6.c -o quintuples6 -O2 -I${HOME}/mpfi-1.5.1/src -I${HOME}/mpfr-3.1.2/src -L${HOME}/mpfi-1.5.1/src/.libs -lmpfi

ineq:		ineq.cpp ../includes/int_double12.0.h 
	g++ ineq.cpp -o ineq -I${HOME}/crlibm-1.0beta4/ -L${HOME}/crlibm-1.0beta4/ -fomit-frame-pointer -frounding-math -finline-functions -mfpmath=387 -lcrlibm

clean:
	rm -f quintuples7 quintuples6 ineq


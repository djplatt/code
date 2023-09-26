IEEE_FLAGS=-O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${HOME}/local/include -I ${HOME}/crlibm-master/ -I ${HOME}/primesieve-5.4.2/include -L ${HOME}/local/lib -L ${HOME}/crlibm-master/ -L ~/primesieve-5.4.2/.libs -lprimesieve -lcrlibm -lm  -lpari
ARB_FLAGS=-O3 -std=c99 -I ${HOME}/local/include -I ${HOME}/primesieve-5.4.2/include/ -I /data/local/include -L ${HOME}/local/lib -L ${HOME}/primesieve-5.4.2/.libs -L /data/local/lib -larb -lprimesieve -lpari
all:	ca1 ca2 ca2_ieee test Rt primorials_ieee

ca1:	ca1.c pari.c
	gcc ca1.c -o ca1 ${ARB_FLAGS} 

ca2:	ca2.c pari.c
	gcc ca2.c -o ca2 ${ARB_FLAGS}

ca_ieee:	ca_ieee.cpp
	g++ ca2_ieee.cpp -o ca2_ieee ${IEEE_FLAGS} 

ca2_ieee:	ca2_ieee.cpp ../includes/int_double14.0.h pari.c
	g++ ca2_ieee.cpp -o ca2_ieee ${IEEE_FLAGS}

primorials_ieee:	primorials_ieee.cpp
	 g++ primorials_ieee.cpp -o primorials_ieee ${IEEE_FLAGS}

Rt:	Rt.cpp ../includes/int_double14.0.h
	g++ Rt.cpp -o Rt ${IEEE_FLAGS}

test:	test.cpp ../includes/int_double14.0.h
	g++ test.cpp -o test ${IEEE_FLAGS}

clean:
	rm -f ca1 ca2 ca_ieee ca2_ieee Rt primorials_ieee

INCLUDE=-I/data/local/include -I${HOME}/primesieve-5.4.2/include -I${HOME}/smalljac_v4.1.3
LIBRARY=-L/data/local/lib -L${HOME}/primesieve-5.4.2/.libs
INCLUDE_G2=-I${HOME}/local/include -I/data/local/include -I${HOME}/primesieve-5.4.2/include
LIBRARY_G2=-L${HOME}/local/lib -L/data/local/lib -L${HOME}/primesieve-5.4.2/.libs
#HWLPOLY=
HWLPOLY=-lhwlpoly

all: fact

fact:	fact.cpp ../andy_elliptic/pari.c
	echo "Building with includes/libs from /data/local"
	g++ -O2 fact.cpp -o fact -finline-functions -fomit-frame-pointer -fopenmp ${INCLUDE} ${LIBRARY} -lgmp -lflint -larb -lm -lsmalljac -lff_poly -lpari -lprimesieve ${HWLPOLY} -lmpfr



clean: 
	rm -f fact


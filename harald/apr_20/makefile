INC=-I ${HOME}/code/includes -I ${HOME}/primesieve-5.4.1/include -I ${HOME}/arb-2.17.0 -I ${HOME}/flint-2.5.2/ -I ${HOME}/mpfr-3.1.2/src -I ${HOME}/gmp-5.1.3
LIB=-L ${HOME}/primesieve-5.4.1/.libs -L ${HOME}/arb-2.17.0

all:	residues read_residues residues-rs

residues:	residues.cpp
	g++ -O2 residues.cpp -o residues ${INC} ${LIB} -larb

residues-rs:	residues-rs.cpp
	g++ -O2 residues-rs.cpp -o residues-rs ${INC} ${LIB} -larb

read_residues:	read_residues.c
	gcc -O2 read_residues.c -o read_residues ${INC} ${LIB} -larb


clean:
	rm -f residues read_residues test

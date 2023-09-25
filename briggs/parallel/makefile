all:	parallel

parallel:	parallel.c
	gcc -O2 -o parallel parallel.c -I ${HOME}/primesieve/include -I ${HOME}/arb-2.19.0 -I ${HOME}/flint-2.6.3/ -L ${HOME}/arb-2.19.0 -larb -L ${HOME}/primesieve -lprimesieve

clean:
	rm -f parallel

all:	ca2

ca2:	ca2.c pari.c
	gcc -O3 -std=c99 ca2.c -o ca2 -I ${HOME}/local/include -I ${HOME}/local/include/pari -L${HOME}/local/lib -larb -lprimesieve -lpari

clean:
	rm -f ca2

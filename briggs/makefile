all:	ca2

ca2:	ca2.c pari.c
	gcc -O3 -std=c99 ca2.c -o ca2 -I ${HOME}/primesieve-5.4.2/include/ -I /data/local/include -L ${HOME}/primesieve-5.4.2/.libs -L /data/local/lib -larb -lprimesieve -lpari

clean:
	rm -f ca2

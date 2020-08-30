ARB_INC=/data/local/include
ARB_LIB=/data/local/lib64
CC=gcc
CFLAGS=-O2 -fPIC -I${ARB_INC} -L${ARB_LIB} -larb -lmpfr
 

all: kacz

kacz: kacz.c
	${CC} kacz.c ../quad/quad.o -o kacz ${CFLAGS}

clean:
	rm -f kacz
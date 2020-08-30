ARB_INC=/data/local/include
ARB_LIB=/data/local/lib64
CC=gcc
CFLAGS=-O2 -c -fPIC -I${ARB_INC}
OBJ=quad.o 

all: lib

lib: $(OBJ)
	${CC} -shared -o libquad.so ${OBJ} -L ${ARB_LIB} -larb

%.o: %.c ${DEPS}
	${CC} ${CFLAGS} $<


clean:
	rm -f *.o
	rm -f libquad.so
ARB_INC="/usr/local/include"
ARB_LIB="/usr/local/lib"
CC=gcc
CFLAGS=-O2 -c -fPIC -I${ARB_INC}
OBJ=quad.o 

all: lib

lib: $(OBJ)
	${CC} -shared -o libquad.so.1.1 -Wl,-soname,libquad.so.1 ${OBJ} -L ${ARB_LIB} -lflint -lmpfr -lgmp

%.o: %.c ${DEPS}
	${CC} ${CFLAGS} $<


clean:
	rm -f *.o
	rm -f libquad.so

CRDIR=${HOME}/crlibm
PSINC=${HOME}/primesieve-5.4.2/include
PSLIB=${HOME}/primesieve-5.4.2/.libs
all: tildem

tildem: tildem.cpp
	g++ tildem.cpp -otildem -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve


clean:
	rm -f tildem


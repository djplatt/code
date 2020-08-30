CRDIR=${HOME}/crlibm
PSINC=${HOME}/primesieve-5.4.2/include
PSLIB=${HOME}/primesieve-5.4.2/.libs
all: Hv

Hv: Hv.cpp
	g++ Hv.cpp -oHv -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve


clean:
	rm -f Hv


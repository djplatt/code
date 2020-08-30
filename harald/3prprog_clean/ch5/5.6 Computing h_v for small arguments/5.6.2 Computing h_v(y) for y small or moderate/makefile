CRDIR=${HOME}/crlibm
PSINC=${HOME}/primesieve-5.4.2/include
PSLIB=${HOME}/primesieve-5.4.2/.libs
all: hvofy hboldvofy

hvofy: hvofy.cpp
	g++ hvofy.cpp -ohvofy -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve

hboldvofy: hboldvofy.cpp
	g++ -O3 hboldvofy.cpp -ohboldvofy -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve


clean:
	rm -f hvofy
	rm -f hboldvofy



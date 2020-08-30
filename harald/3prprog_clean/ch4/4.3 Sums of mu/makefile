CRDIR=${HOME}/crlibm
PSINC=${HOME}/primesieve-5.4.2/include
PSLIB=${HOME}/primesieve-5.4.2/.libs

all: checkm checkmopt musumE

checkm: checkm.cpp
	g++ checkm.cpp -ocheckm -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve

checkmopt: checkmopt.cpp
	g++ checkmopt.cpp -ocheckmopt -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve

musumE: musumE.cpp
	g++ musumE.cpp -omusemE -fopenmp -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve


clean:
	rm -f checkm
	rm -f checkmopt
	rm -f musumE

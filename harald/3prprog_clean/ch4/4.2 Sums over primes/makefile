CRDIR=${HOME}/crlibm
PSINC=${HOME}/primesieve-5.4.2/include
PSLIB=${HOME}/primesieve-5.4.2/.libs
all: estpsi estpsin esttheta esttheta2 checkm

estpsi: estpsi.cpp
	g++ estpsi.cpp -oestpsi -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve

estpsin: estpsin.cpp
	g++ estpsin.cpp -oestpsin -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve

esttheta: esttheta.cpp
	g++ esttheta.cpp -oesttheta -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve

esttheta2: esttheta2.cpp
	g++ esttheta2.cpp -oesttheta2 -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve

checkm: checkm.cpp
	g++ checkm.cpp -ocheckm -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ${PSINC} -L${PSLIB} -I${CRDIR} -L${CRDIR} -lcrlibm -lprimesieve

clean:
	rm -f estpsi
	rm -f estpsin

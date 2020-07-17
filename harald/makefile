CRDIR=/panfs/panasas01/math/madjp/crlibm-1.0beta4
PSDIR=/panfs/panasas01/math/madjp/primesieve-5.4.1/include
PSLIB=/panfs/panasas01/math/madjp/primesieve-5.4.1/.libs

musumv: musumv7.cpp
	g++ -fopenmp musumv7.cpp -omusumv -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I${CRDIR} -L${CRDIR} -lcrlibm -lgomp

musum: djp_musim.cpp
	g++ djp_musim.cpp -odjp_musim -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I${CRDIR} -I${PSDIR} -L${CRDIR} -lcrlibm -L${PSLIB} -lprimesieve

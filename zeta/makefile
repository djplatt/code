INCLUDE=/data/local/include
LIBRARY=/data/local/lib
all: arb_zeta

arb_zeta: arb_zeta.cpp arb_fft.h
	g++ -O1 arb_zeta.cpp -o arb_zeta -finline-functions -fomit-frame-pointer -I ${INCLUDE} -L ${LIBRARY} -larb -lflint -lgmp


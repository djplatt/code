all: arb_zeta

arb_zeta: arb_zeta.cpp parameters.h
	g++ arb_zeta.cpp -o arb_zeta -O2 -finline-functions -fomit-frame-pointer -I ${HOME}/local/include -L ${HOME}/local/lib -lflint

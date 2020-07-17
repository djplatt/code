sieve: sieve.cpp
	g++ -fPIC -o sieve sieve.cpp -I ${HOME}/primesieve-5.4.1/include -O2 -finline-functions -fomit-frame-pointer -L ${HOME}/primesieve-5.4.1/.libs -lprimesieve

theta: theta.cpp ../includes/int_double12.0.h
	g++ -fPIC -o theta theta.cpp -I${HOME}/crlibm-1.0beta4/ -L${HOME}/crlibm-1.0beta4/ -fomit-frame-pointer -frounding-math -finline-functions -mfpmath=387 -lcrlibm

theta2: theta2.cpp ../includes/int_double12.0.h
	g++ -o theta2 theta2.cpp -I${HOME}/crlibm-1.0beta4/ -I ${HOME}/primesieve-5.4.1/include -L${HOME}/crlibm-1.0beta4/ -fomit-frame-pointer -frounding-math -finline-functions -mfpmath=387 -lcrlibm -L ${HOME}/primesieve-5.4.1/.libs -lprimesieve

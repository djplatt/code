all: siegel

siegel: siegel.cpp ../includes/int_double12.0.h ../characters-master/slint.h 
	g++ -o siegel siegel.cpp -I${HOME}/crlibm-1.0beta4/ -L${HOME}/crlibm-1.0beta4/ -fomit-frame-pointer -frounding-math -finline-functions -mfpmath=387 -lcrlibm

clean:
	rm -f siegel


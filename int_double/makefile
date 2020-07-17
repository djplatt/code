all: test_int_double

test_int_double: test_int_double.cpp int_double.h
	g++ -O3 -frounding-math -finline-functions -fomit-frame-pointer -o test_int_double test_int_double.cpp -I ${HOME}/crlibm-1.0beta4/ -L ${HOME}/crlibm-1.0beta4/ -lcrlibm


clean:
	rm -f test_int_double
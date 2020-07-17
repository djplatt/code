all: test_log1p test_abs

test_log1p: test_log1p.cpp ../includes/int_double13.0.h
	g++ test_log1p.cpp -o test_log1p -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm

test_abs: test_abs.cpp ../includes/int_double14.0.h
	g++ test_abs.cpp -o test_abs -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I ~/crlibm-1.0beta4/ -L ~/crlibm-1.0beta4/ -lcrlibm -lm

clean:
	rm -f test_log1p


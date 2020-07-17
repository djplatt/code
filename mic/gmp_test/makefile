all: test

test: test.c
	/opt/intel/bin/icc -mmic -o test test.c -I${HOME}/gmp-5.1.3-mic/ -L${HOME}/gmp-5.1.3-mic/.libs/ -lgmp
	echo "now sftp <test> to host mic0"
	echo "then ssh mic0 ./test"

clean:
	rm -f test



all: hello

hello: hello.c
	icc -mmic -o hello hello.c
	echo "now sftp <hello> to host mic0"
	echo "then ssh mic0 ./hello"

clean:
	rm -f hello



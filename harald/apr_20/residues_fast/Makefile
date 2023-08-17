all: residues

residues: residues.cpp parameters.h
	g++ residues.cpp -o residues -O2 -finline-functions -fomit-frame-pointer -lflint -larb

clean:
	rm -f residues

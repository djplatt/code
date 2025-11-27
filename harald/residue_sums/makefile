INCDIR="."
LIBDIR="."
all:	residue_sums_nofft arb_load_sum

residue_sums_nofft:	residue_sums_nofft.c
	gcc -O2 -o residue_sums_nofft -I${INCDIR} -L ${LIBDIR} residue_sums_nofft.c -lflint

arb_load_sum:	arb_load_sum.c
	gcc -O2 -o arb_load_sum -I${INCDIR} -L ${LIBDIR} arb_load_sum.c -lflint 

clean:
	rm -rf residue_sums_nofft arb_load_sum

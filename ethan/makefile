all:	check_first_zero sum_gamma sum_rho sum_jun_23

check_first_zero:	check_first_zero.c
	gcc -O2 -std=c99 check_first_zero.c -o check_first_zero -I ${HOME}/arb-2.19.0/ -I ${HOME}/flint-2.6.3 -L ${HOME}/arb-2.19.0/ -larb

sum_rho:	sum_rho.c
	gcc -O2 -o sum_rho sum_rho.c -I ${HOME}/arb-2.19.0/ -I ${HOME}/flint-2.6.3 -L ${HOME}/arb-2.19.0/ -larb

sum_gamma:	sum_gamma.c
	gcc -O2 -o sum_gamma sum_gamma.c -I ${HOME}/arb-2.19.0/ -I ${HOME}/flint-2.6.3 -L ${HOME}/arb-2.19.0/ -larb

sum_jun_23:	sum_jun_23.c
	gcc -O2 -o sum_jun_23 sum_jun_23.c -I ${HOME}/arb-2.19.0/ -I ${HOME}/flint-2.6.3 -L ${HOME}/arb-2.19.0/ -larb

clean:
	rm check_first_zero sum_rho sum_gamma sum_jun_23

#ifndef PARAMETERS
#define PARAMETERS
//
// MATHEMATICS OF COMPUTATION
// Volume 86, Number 307, September 2017, Pages 2449–2467
// http://dx.doi.org/10.1090/mcom/3198
// Article electronically published on February 13, 2017



// Parameters for arb_zeta

// limits for T0
#define T0_MAX ((double) 3.0100000000000000000000000000000000000e12)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)

// we do a cheeky upsample of this factor
#define UPSAM ((int) 4)

// how many points. Power of 2 please
#define N ((int) 131072) // steps 6-8
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2) // f/F spacing
#define N1 ((int) N/UPSAM) // intermediate FFT length steps 1-6
#define one_over_A1 ((double) one_over_A*UPSAM) // g/G spacing
#define B ((double) N*one_over_A) // 5376

#define h ((double) 1.1600000000000000000000000000000000000e2) // Gaussian decay in (3.1)

// lemma A.11 
#define ftwiderr_d ((double) 2.37e-115)

// lemma A.9
#define fhattwiderr_d ((double) 1e-307)

#define M ((int) 768000) // lemma B.1 (called J)
#define fhatsumerr_d ((double) 3.26e-33)

#define K ((int) 23) // lemma B.2
#define tayerr_d ((double) 6.73e-33) // lemma B.2

// lemma A.5
#define gtwiderr_d ((double) 1.83e-44)

//lemma A.7
#define Gtwiderr_d ((double) 1e-307)

// lemma A.8
#define Fmaxerr_d ((double) 3.93e-245)

// upsampling better than this if one uses O(t^(-1/4)) estimate
// for gamma(1/4+it/2)exp(Pi t/4)
//
// Upsampling parameters
#define INTER_SPACING ((int) 1) // allows for skipping points. not used.
#define INTER_A ((double) INTER_SPACING*one_over_A)
#define H ((double) 2.0312500000000000000000000000000000000e-1) // in defn of W(t)
#define Ns ((int) 70) // lemma C.3
#define intererr_d ((double) 2.45e-40) // lemma C.3
//************************************

#define EXTRA_PREC (50) // more than log2(t0/(2Pi)*log(Msqrt(Pi)))

// Turing method parameters
#define TURING_WIDTH ((double) 21.0)
#define TURING_LEN ((int) (TURING_WIDTH/one_over_A))

// Newton Iteration Parameters
#define NEWTON_ITS (8) // Do this many iterations with Newton Raphson to locate roots

#endif


// see zeta-fft-errors.gp
// parameters for t0<=3.062*10^10
#define T0_MAX (3.062e10) // maximum t0
#define EXTRA_PREC (50) // more than log2(t0/(2Pi)*log(Msqrt(Pi)))
#define T0_MIN (5000.0) // mimimal input t0
#define UPSAM (32) // the rate of upsampling from g to f
#define N ((int) 1<<20) // final FFT length
#define N1 ((int) N/UPSAM) // intermediate FFT length
#define one_over_A ((double) 21.0/4096.0) // output spacing
#define one_over_A1 ((double) one_over_A*UPSAM) // g/G spacing
#define h ((double) 176431.0/2048.0) // Gaussian width for f(t)
#define B ((double) N*one_over_A)
#define M ((int) 104000) // number of terms in F(x) sum

#define K ((int) 44) // number of Taylor terms
#define gtwiderr_d ((double) 3.2e-82)
#define Gtwiderr_d ((double) 6.1e-213)
//#define fhatsumerr_d ((double) 4.3e-83)
#define fhatsumerr_d ((double) 1.28e-83) // New version sig=1675
#define tayerr_d ((double) 1.5e-82)
#define fhattwiderr_d ((double) 1.0e-307) // actually much smaller
#define ftwiderr_d ((double) 8.1e-211)
#define Fmaxerr_d ((double) 1.0e-307) // Max mod of F(N1/(2B)) actually much smaller

// Turing method parameters
#define TURING_WIDTH (42)
#define TURING_LEN ((int) (TURING_WIDTH/one_over_A))

// upsampling better than this if one uses O(t^(-1/4)) estimate
// for gamma(1/4+it/2)exp(Pi t/4)
//
//#define UPSAM_HIGH
#ifndef UPSAM_HIGH
// Upsampling parameters
#define Ns (70) // take this many points either side of t0
#define H ((double) 2089.0/16384.0) // Gaussian width for upsampling
#define INTER_SPACING (5) // Use every 5th sample
#define intererr_d ((double) 5.0e-41)
#define dintererr_d ((double) 1e-30)
#define INTER_A ((double) INTER_SPACING*one_over_A)
#else
#define Ns (140) // take this many points either side of t0
#define H ((double) 23855.0/(1<<17)) // Gaussian width for upsampling
#define INTER_SPACING (5) // Use every 5th sample
#define intererr_d ((double) 3.7e-83)
#define INTER_A ((double) INTER_SPACING*one_over_A)
#endif

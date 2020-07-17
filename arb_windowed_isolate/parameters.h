/*
// parameters for t0<=3*10^10
// see zeta-fft-errors.gp
#define T0_MAX (3.0e10) // maximum t0
#define EXTRA_PREC (50) // more than log2(t0/(2Pi)*log(Msqrt(Pi)))
#define T0_MIN (5000.0) // mimimal input t0
#define UPSAM (32) // the rate of upsampling from g to f
#define N ((int) 1<<20) // final FFT length
#define N1 ((int) N/UPSAM) // intermediate FFT length
#define one_over_A ((double) 21.0/4096.0) // output spacing
#define one_over_A1 ((double) one_over_A*UPSAM) // g/G spacing
#define h ((double) 176431.0/2048.0) // Gaussian width for f(t)
#define B ((double) N*one_over_A)
#define M ((int) 103000) // number of terms in F(x) sum

#define K ((int) 44) // number of Taylor terms
#define gtwiderr_d ((double) 3.2e-82)
#define Gtwiderr_d ((double) 2.3e-213)
#define fhatsumerr_d ((double) 1.7e-83)
#define tayerr_d ((double) 1.5e-82)
#define fhattwiderr_d ((double) 1.0e-307) // actually much smaller
#define ftwiderr_d ((double) 8.1e-211)
#define Fmaxerr_d ((double) 1.0e-307) // Max mod of F(N1/(2B)) actually much smaller

// Turing method parameters
#define TURING_WIDTH (42)
#define TURING_LEN ((int) (TURING_WIDTH/one_over_A))

// Upsampling parameters
#define Ns (70) // take this many points either side of t0
#define H ((double) 2089.0/16384.0) // Gaussian width for upsampling
#define INTER_SPACING (5) // Use every 5th sample
#define intererr_d ((double) 5.0e-41)
#define INTER_A ((double) INTER_SPACING*one_over_A)

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
#define INTER_A ((double) INTER_SPACING*one_over_A)
#else
#define Ns (140) // take this many points either side of t0
#define H ((double) 23855.0/(1<<17)) // Gaussian width for upsampling
#define INTER_SPACING (5) // Use every 5th sample
#define intererr_d ((double) 3.7e-83)
#define INTER_A ((double) INTER_SPACING*one_over_A)
#endif
*/

// test parameters for k=23
/*
#define T0_MAX (1.001e13) // maximum t0
#define EXTRA_PREC (50) // more than log2(t0/(2Pi)*log(Msqrt(Pi)))
#define T0_MIN (1.0e10) // mimimal input t0
#define UPSAM (32) // the rate of upsampling from g to f
#define N ((int) 1<<20) // final FFT length
#define N1 ((int) N/UPSAM) // intermediate FFT length
#define one_over_A ((double) 21.0/4096.0) // output spacing
#define one_over_A1 ((double) one_over_A*UPSAM) // g/G spacing
#define h ((double) 116) // Gaussian width for f(t)
#define B ((double) N*one_over_A)
#define M ((int) 1670000) // number of terms in F(x) sum

#define K ((int) 23) // number of Taylor terms
#define gtwiderr_d ((double) 1.9e-44)
#define Gtwiderr_d ((double) 1e-307)
#define fhatsumerr_d ((double) 1.1e-34)
#define tayerr_d ((double) 1.0e-32)
#define fhattwiderr_d ((double) 1e-307) // actually smaller 
#define ftwiderr_d ((double) 3.0e-115)
#define Fmaxerr_d ((double) 1.1e-211) // Max mod of F(N1/(2B)) 

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
#define intererr_d ((double) 1.31e-38)
#define INTER_A ((double) INTER_SPACING*one_over_A)
#else

#define Ns (140) // take this many points either side of t0
#define H ((double) 23855.0/(1LL<<17)) // Gaussian width for upsampling
#define INTER_SPACING (5) // Use every 5th sample
#define intererr_d ((double) 5.76e-81)
#define INTER_A ((double) INTER_SPACING*one_over_A)
#endif
*/
// test parameters for k=24
/*
#define T0_MAX (1.001e13) // maximum t0
#define EXTRA_PREC (50) // more than log2(t0/(2Pi)*log(Msqrt(Pi)))
#define T0_MIN (1.0e10) // mimimal input t0
#define UPSAM (32) // the rate of upsampling from g to f
#define N ((int) 1<<20) // final FFT length
#define N1 ((int) N/UPSAM) // intermediate FFT length
#define one_over_A ((double) 21.0/4096.0) // output spacing
#define one_over_A1 ((double) one_over_A*UPSAM) // g/G spacing
#define h ((double) 114) // Gaussian width for f(t)
#define B ((double) N*one_over_A)
#define M ((int) 1680000) // number of terms in F(x) sum

#define K ((int) 18)//24) // number of Taylor terms
#define gtwiderr_d ((double) 9.99e-46)
#define Gtwiderr_d ((double) 1e-307)
#define fhatsumerr_d ((double) 2.9e-37)
#define tayerr_d ((double) 8.82e-35)
#define fhattwiderr_d ((double) 1e-307) // actually smaller 
#define ftwiderr_d ((double) 2.2e-119)
#define Fmaxerr_d ((double) 1.1e-211) // Max mod of F(N1/(2B)) 

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
#define intererr_d ((double) 1.31e-38)
#define INTER_A ((double) INTER_SPACING*one_over_A)
#else

#define Ns (140) // take this many points either side of t0
#define H ((double) 23855.0/(1LL<<17)) // Gaussian width for upsampling
#define INTER_SPACING (5) // Use every 5th sample
#define intererr_d ((double) 5.76e-81)
#define INTER_A ((double) INTER_SPACING*one_over_A)

#endif


// Newton Iteration Parameters
#define NEWTON_ITS (8) // Do this many iterations with Newton Raphson to locate roots



#define T0_MAX (1.001e13) // maximum t0
#define EXTRA_PREC (50) // more than log2(t0/(2Pi)*log(Msqrt(Pi)))
#define T0_MIN (1.0e10) // mimimal input t0
#define UPSAM (32) // the rate of upsampling from g to f
#define N ((int) 1<<20) // final FFT length
#define N1 ((int) N/UPSAM) // intermediate FFT length
#define one_over_A ((double) 21.0/4096.0) // output spacing
#define one_over_A1 ((double) one_over_A*UPSAM) // g/G spacing
#define h ((double) 140) // Gaussian width for f(t)
#define B ((double) N*one_over_A)
#define M ((int) 1650000) // number of terms in F(x) sum

#define K ((int) 18) // number of Taylor terms
#define gtwiderr_d ((double) 1.5e-22)
#define Gtwiderr_d ((double) 1e-317) // smaller in truth
#define fhatsumerr_d ((double) 2.9e-39)
#define tayerr_d ((double) 5.3e-22)
#define fhattwiderr_d ((double) 1e-307) // actually smaller 
#define ftwiderr_d ((double) 1.6e-78)
#define Fmaxerr_d ((double) 1.2e-211) // Max mod of F(N1/(2B)) 

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
#define intererr_d ((double) 1.31e-38)
#define INTER_A ((double) INTER_SPACING*one_over_A)
#else

#define Ns (140) // take this many points either side of t0
#define H ((double) 23855.0/(1LL<<17)) // Gaussian width for upsampling
#define INTER_SPACING (5) // Use every 5th sample
#define intererr_d ((double) 5.76e-81)
#define INTER_A ((double) INTER_SPACING*one_over_A)

#endif


// Newton Iteration Parameters
#define NEWTON_ITS (8) // Do this many iterations with Newton Raphson to locate roots
*/

/*
//main(21/4096,2^20,32,1e11,23,5,70,2089/16384,141000,116)
#define T0_MAX ((double) 10.000000000000000000000000000000000000e10)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 32)
#define N ((int) 1048576)
#define one_over_A ((double) 5.1269531250000000000000000000000000000e-3)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.20e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 141000)
#define fhatsumerr_d ((double) 8.62e-38)
#define K ((int) 23)
#define tayerr_d ((double) 7.20e-36)
#define gtwiderr_d ((double) 4.57e-47)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 1e-307)
#define INTER_SPACING ((int) 5)
#define H ((double) 1.2750244140625000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 3.55e-41)
//************************************
*/

/*
//main(21/2048,2^19,16,1e11,23,2,70,2089/16384,141000,116)
#define T0_MAX ((double) 10.000000000000000000000000000000000000e10)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 16)
#define N ((int) 524288)
#define one_over_A ((double) 1.0253906250000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.20e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 141000)
#define fhatsumerr_d ((double) 8.62e-38)
#define K ((int) 23)
#define tayerr_d ((double) 7.20e-36)
#define gtwiderr_d ((double) 4.57e-47)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 1e-307)
#define INTER_SPACING ((int) 2)
#define H ((double) 1.2750244140625000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 4.07e-26)
//************************************
*/

/*
//main(21/1024,2^18,8,1e11,23,1,70,2089/16384,141000,116)
#define T0_MAX ((double) 10.000000000000000000000000000000000000e10)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 8)
#define N ((int) 262144)
#define one_over_A ((double) 2.0507812500000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.20e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 141000)
#define fhatsumerr_d ((double) 8.62e-38)
#define K ((int) 23)
#define tayerr_d ((double) 7.20e-36)
#define gtwiderr_d ((double) 4.57e-47)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 1e-307)
#define INTER_SPACING ((int) 1)
#define H ((double) 1.2750244140625000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 4.07e-26)
//************************************
*/

/*
//main(21/512,2^17,4,1e11,23,1,70,13/64,141000,116)
#define T0_MAX ((double) 10.000000000000000000000000000000000000e10)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.20e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 141000)
#define fhatsumerr_d ((double) 8.62e-38)
#define K ((int) 23)
#define tayerr_d ((double) 2.89e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 1e-307)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 1.09e-41)
//************************************
*/

/*
//main(21/512,2^17,4,5e10,23,1,70,13/64,99000,116)
#define T0_MAX ((double) 5.0000000000000000000000000000000000000e10)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.05e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 99000)
#define fhatsumerr_d ((double) 2.58e-33)
#define K ((int) 23)
#define tayerr_d ((double) 2.42e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 1e-307)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 6.77e-42)
//************************************
*/

/*
//main(21/256,2^16,2,1e11,23,1,85,1/2,141000,116)
#define T0_MAX ((double) 10.000000000000000000000000000000000000e10)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 2)
#define N ((int) 65536)
#define one_over_A ((double) 8.2031250000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.20e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 141000)
#define fhatsumerr_d ((double) 8.62e-38)
#define K ((int) 23)
#define tayerr_d ((double) 2.89e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 1e-307)
#define INTER_SPACING ((int) 1)
#define H ((double) 5.0000000000000000000000000000000000000e-1)
#define Ns ((int) 85)
#define intererr_d ((double) 1.83e-41)
//************************************
*/
/*
//main(21/512,2^17,4,2e11,23,1,70,13/64,198000,116)
#define T0_MAX ((double) 2.0000000000000000000000000000000000000e11)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.38e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 198000)
#define fhatsumerr_d ((double) 2.58e-33)
#define K ((int) 23)
#define tayerr_d ((double) 3.42e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 7.84e-303)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 1.90e-41)
//************************************
*/

/*
//main(21/512,2^17,4,3e11,23,1,70,13/64,250000,116)
#define T0_MAX ((double) 3.0000000000000000000000000000000000000e11)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.50e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 250000)
#define fhatsumerr_d ((double) 1.13e-54)
#define K ((int) 23)
#define tayerr_d ((double) 3.84e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 3.34e-294)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 2.71e-41)
//************************************
*/
/*
//main(21/512,2^17,4,5e11,23,1,70,13/64,313000,116) 
#define T0_MAX ((double) 5.0000000000000000000000000000000000000e11)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.66e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 313000)
#define fhatsumerr_d ((double) 3.46e-33)
#define K ((int) 23)
#define tayerr_d ((double) 4.30e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 2.48e-283)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 4.32e-41)
//************************************
*/

/*
// main(21/512,2^17,4,6e11,23,1,70,13/64,343000,116) 
#define T0_MAX ((double) 6.0000000000000000000000000000000000000e11)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.72e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 343000)
#define fhatsumerr_d ((double) 2.07e-33)
#define K ((int) 23)
#define tayerr_d ((double) 4.50e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 1.88e-279)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 5.12e-41)
//************************************
*/
/*
//main(21/512,2^17,4,7e11,23,1,70,13/64,371000,116)
#define T0_MAX ((double) 7.0000000000000000000000000000000000000e11)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.78e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 371000)
#define fhatsumerr_d ((double) 2.84e-34)
#define K ((int) 23)
#define tayerr_d ((double) 4.68e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 3.58e-276)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 5.92e-41)
//************************************
*/
/*
//main(21/512,2^17,4,8e11,23,1,70,13/64,396000,116)
#define T0_MAX ((double) 8.0000000000000000000000000000000000000e11)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.82e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 396000)
#define fhatsumerr_d ((double) 2.58e-33)
#define K ((int) 23)
#define tayerr_d ((double) 4.83e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 2.49e-273)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 6.72e-41)
//************************************


//main(21/512,2^17,4,9e11,23,1,70,13/64,420000,116)
#define T0_MAX ((double) 9.0000000000000000000000000000000000000e11)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.87e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 420000)
#define fhatsumerr_d ((double) 2.77e-33)
#define K ((int) 23)
#define tayerr_d ((double) 4.98e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 7.98e-271)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 7.52e-41)
//************************************
*/
/*
//main(21/512,2^17,4,1.2e12,23,1,70,13/64,485000,116) 
#define T0_MAX ((double) 1.2000000000000000000000000000000000000e12)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 1.98e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 485000)
#define fhatsumerr_d ((double) 2.58e-33)
#define K ((int) 23)
#define tayerr_d ((double) 5.35e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 1.07e-264)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 9.93e-41)
//************************************


//main(21/512,2^17,4,1.4e12,23,1,70,13/64,524000,116) 
#define T0_MAX ((double) 1.4000000000000000000000000000000000000e12)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 2.04e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 524000)
#define fhatsumerr_d ((double) 1.77e-33)
#define K ((int) 23)
#define tayerr_d ((double) 5.56e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 2.02e-261)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 1.16e-40)
//************************************
*/
//main(21/512,2^17,4,1.6e12,23,1,70,13/64, 560000,116)
#define T0_MAX ((double) 1.6000000000000000000000000000000000000e12)
#define T0_MIN ((double) 1.0000000000000000000000000000000000000e10)
#define UPSAM ((int) 4)
#define N ((int) 131072)
#define one_over_A ((double) 4.1015625000000000000000000000000000000e-2)
#define h ((double) 1.1600000000000000000000000000000000000e2)
#define ftwiderr_d ((double) 2.09e-115)
#define fhattwiderr_d ((double) 1e-307)
#define M ((int) 560000)
#define fhatsumerr_d ((double) 2.77e-33)
#define K ((int) 23)
#define tayerr_d ((double) 5.75e-33)
#define gtwiderr_d ((double) 1.83e-44)
#define Gtwiderr_d ((double) 1e-307)
#define Fmaxerr_d ((double) 1.41e-258)
#define INTER_SPACING ((int) 1)
#define H ((double) 2.0312500000000000000000000000000000000e-1)
#define Ns ((int) 70)
#define intererr_d ((double) 1.32e-40)
//************************************

#define EXTRA_PREC (50) // more than log2(t0/(2Pi)*log(Msqrt(Pi)))
#define N1 ((int) N/UPSAM) // intermediate FFT length
#define one_over_A1 ((double) one_over_A*UPSAM) // g/G spacing
#define B ((double) N*one_over_A)
// Turing method parameters
#define TURING_WIDTH ((double) 21.0)
#define TURING_LEN ((int) (TURING_WIDTH/one_over_A))

// upsampling better than this if one uses O(t^(-1/4)) estimate
// for gamma(1/4+it/2)exp(Pi t/4)
//
// Upsampling parameters
#define INTER_A ((double) INTER_SPACING*one_over_A)

// Newton Iteration Parameters
#define NEWTON_ITS (8) // Do this many iterations with Newton Raphson to locate roots


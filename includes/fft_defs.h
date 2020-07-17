#ifndef FFT_DEFS
#define FFT_DEFS

#define MAX_Q (400000)
#define MAX_CONV (1<<19)//23) // smallest pwr of 2 >=phi(Q)-1
// this is sufficient because for prime q we do two half length Bluesteins
#define MAX_FACS (7) // 2*3*5*7*11*13*17>MAX_Q
#define MAX_DIMS (MAX_FACS+1) // plus 1 for 2^n trick
#define MAX_SIMPLE_DFT (31) // use simple O(n^2) algorithm for n<= this

// structure to hold factor information
typedef struct
{
	long unsigned int pr;   // = 0 if no primitive root
	long unsigned int phi;  // Euler phi(q)
	long unsigned int num_facs;  // number of prime power factors
	long unsigned int facs[MAX_FACS];  // the factors p^n
	long unsigned int primes[MAX_FACS]; // the prime factors p
} factor;

#endif

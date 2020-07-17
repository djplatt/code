#ifndef FFT_DEFS
#define FFT_DEFS


#define MAX_Q (2000000) //
#define MAX_CONV (1<<22) // smallest pwr of 2 >=2Q-1
#define MAX_FACS (8) // was 9 2*3*5*7*11*13*17*19>MAX_Q
#define MAX_DIMS (MAX_FACS+1) // plus 1 for 2^n trick
#define MAX_SIMPLE_DFT (31) // use simple O(n^2) algorithm for n<= this

// structure to hold factor information
typedef struct
{
	unsigned long int pr;   // = 0 if no primitive root
	unsigned long int phi;  // Euler phi(q)
	unsigned long int num_facs;  // number of prime power factors
	unsigned long int facs[MAX_FACS];  // the factors p^n
	unsigned long int primes[MAX_FACS]; // the prime factors p
} factor;


inline int gcd (unsigned int a, unsigned int b)
/* Euclid algorithm gcd */
{
	unsigned int c;
	while(a!=0)
	{
		c=a;
		a=b%a;
		b=c;
	};
	return(b);
};

inline int co_prime(unsigned int a, unsigned int b)
{
	return(gcd(a,b)==1);
};


#endif

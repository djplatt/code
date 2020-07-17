#ifndef IM_S
#define IM_S
// im_s.h
// defines structure to hold Im(s) and Gamma Values

// structure to hold gamma values etc.
#if defined(MS)
__declspec(align(16)) struct Im_s

{
	int_complex lambda_s; // gamma(s/2)
	int_complex lambda_s_a; // gamma((s+1)/2)
	double im_s; // imaginary part of s
};

typedef Im_s im_s;
#else
typedef struct

{
	int_complex lambda_s; // gamma(s/2)
	int_complex lambda_s_a; // gamma((s+1)/2)
	int_complex pi_minus_it_2; // pi^(-it/2)
	double im_s; // imaginary part of s
} 
__attribute__ ((aligned(16)))
im_s;
#endif
#endif

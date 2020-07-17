#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include <complex>
#include "../includes/int_double8.0.h"
#include "../includes/im_s.h"

#ifndef LINUX
#include <windows.h>
#endif

#define VERSION "2.0"
// 14/1/09
// 1.0 Original
// 2.0 gam_s values now include Pi factor
//     hur values now include q^s/2 factor

/* Andy Booker's Interval FFT code */

void initfft(unsigned int n, int_complex *w) // n a power of 2, w of size n/2
{
	unsigned int k;
	int_double two_pi_n;

	w[0]=c_one;
	two_pi_n=(d_pi*2)/n;
	
	for(k=1;k<(n/2);k++)
		sin_cos(two_pi_n*k,&w[k].imag,&w[k].real);
}

/* perform an in place FFT ,n a power of 2 */ 
void fft(int_complex *x,unsigned int n,int_complex *w)
{
	unsigned int i,j,k,l;
	int_complex *p,*xend=x+n,ctemp;

	/* swap each element with one with bit-reversed index */
	for (i=0,l=n>>1;i<l;++i) {
		/* j = bit reversal of i */
		for (k=1,j=0;k<n;k<<=1) {
			j <<= 1;
			if (i & k)
				j |= 1;
		}
		if (i < j)
		{
			ctemp=x[i];
			x[i]=x[j];
			x[j]=ctemp;
		}
		else if (i > j)
		{
			ctemp=x[n-1-i];
			x[n-1-i]=x[n-1-j];
			x[n-1-j]=ctemp;
		}
		++i, j |= l;
		ctemp=x[i];
		x[i]=x[j];
		x[j]=ctemp;
	}

	for (k=1,l=n/2;k<n;k<<=1,l>>=1)
		for (p=x;p<xend;p+=k)
			for (j=0;j<n/2;j+=l,p++) {
				ctemp=p[k]*w[j];
				p[k]=p[0]-ctemp;
				p[0]+=ctemp;
			}
}

/* perform an in place inverse FFT */
void ifft(int_complex *x,unsigned int n,int_complex *w) {
	unsigned int i,l=n>>1;
	int_complex ctemp;

	fft(x,n,w);
	x[0]=x[0]/n;
	x[l]=x[l]/n;

	for (i=1;i<l;i++) {
		x[i]=x[i]/n;
		x[n-i]=x[n-i]/n;
		ctemp=x[i];
		x[i]=x[n-i];
		x[n-i]=ctemp;
	}
}

/* circular convolution; f and g must be distinct */
void convolve(int_complex *result,int_complex *f,int_complex *g,unsigned int n, int_complex *w)
{
	unsigned int i;
	fft(f,n,w);
	fft(g,n,w);
	for (i=0;i<n;i++)
		result[i]=f[i]*g[i];
	ifft(result,n,w);
}

#define MAX_N (100000)  // largest to 100,000 will be 99,990
#define TWO_MAX_N (MAX_N*2)

int_complex w[MAX_N],a[TWO_MAX_N],b[TWO_MAX_N],b_spare[TWO_MAX_N],b_star_conj[TWO_MAX_N];

void init_bluestein_fft(unsigned int n, unsigned int *conv_size,int_complex *w,
					   int_complex *b, int_complex *b_star_conj)
/*
   n is an unsigned integer <=MAX_N
   conv_size[0] will contain smallest power of 2 >= 2n-1
   w will contain conv_size/2 powers of exp(2*pi*i/conv_size)
   b will contain conv_size (zero padded) values of exp(pi*i*k^2/n)
   b_star_conj will contain n values of conj(b[k])
   */
{
	unsigned int i,isqrmod;
	unsigned long long ismod;
	int_double pi_over_n;

	/* lets make sure it isnt
	if(n>MAX_N)
	{
		printf("Fatal error, n too large (%d) in init_bluestein_fft. Exiting.\n",n);
		exit(0);
	};
	*/

	conv_size[0]=1;
	while(conv_size[0]<2*n-1)
		conv_size[0]<<=1;
	//printf("FFT size set to %d.\n",conv_size);
	initfft(conv_size[0],w);
	pi_over_n=d_pi/n;
	b[0]=c_one;
	b_star_conj[0]=c_one;
	for(i=1;i<n;i++)
	{
		ismod=(long long unsigned int)i*(long long unsigned int)i;
		ismod=ismod%(2*n);
		isqrmod=(unsigned int) ismod;
		
		sin_cos(pi_over_n*isqrmod,&b[i].imag,&b[i].real);
		
		b_star_conj[i]=conj(b[i]);
	}
	/*
	int r_e,worst=100;
	for(i=1;i<n;i++)
	{
		r_e=-rel_error(b[i].real);
		if(r_e<worst)
			worst=r_e;
	}
	printf("%d %d\n",n,worst);
	*/

	for(i=1;i<n;i++)
		b[conv_size[0]-i]=b[i];
	for(i=n;i<=(conv_size[0]-n);i++)
		b[i]=c_zero;

/*
	for(i=1;i<n;i<<=1)
	{
		printf("w^(%d^2)=",i);
		print_int_complex(b[i]);
		printf("\n");
	}
*/
}



void bluestein_fft(unsigned int stride, int_complex *x, unsigned int n,
				   int_complex *w, int_complex *a, int_complex *b,
				   int_complex *b_star_conj, unsigned int conv_size,
				   int_complex *b_spare)
/* 
   x contains n inputs spaced out by stride starting at x[0]
   w[k] contains exp(2*pi*i*k/conv_size) for k=0..(conv_size/2)-1
   a[k] contains rubbish for k=0..conv_size-1
   b[k] contains exp(pi*i*k^2/n) for k=0..n-1
   b[conv_size-k] contains b[k] for k=1..n-1
   b[k] contains zero otherwise
   b_star_conj[k] contains conj(b[k]) for k=0..n-1
   conv_size is smallest power of 2 >=2n-1
   */
{
	unsigned int i,j;

	a[0]=x[0]; 
	for(i=1,j=stride;i<n;i++,j+=stride)
		a[i]=x[j]*b_star_conj[i];
	for(i=n;i<conv_size;i++)
		a[i]=c_zero;
	for(i=0;i<conv_size;i++)
		b_spare[i]=b[i];
	convolve(a,a,b_spare,conv_size,w); // a<- a*b
	x[0]=a[0];
	for(i=1,j=stride;i<n;i++,j+=stride)
		x[j]=a[i]*b_star_conj[i];
}

#define MAX_FACS (6)  // 2*4*3*5*7*11*13 > 100,000
// FFT_VEC_SIZE must be at least as large as max(phi(q))
#define FFT_VEC_SIZE (512*500)  // each one is 32 bytes for double complex 

int_complex fft_vec[FFT_VEC_SIZE];

typedef int_complex dcomplex;

// structure to hold factor information
typedef struct
{
  unsigned int pr;   // = 0 if no primitive root
  unsigned int phi;  // Euler phi(q)
  unsigned int num_facs;  // number of prime power factors
  unsigned int facs[MAX_FACS];  // the factors p^n
  unsigned int primes[MAX_FACS]; // the prime factors p
} factor;


/*
// Complex Gamma Function
// see http://bytes.com/forum/thread576697.html 
static const int g=7;
static const double pi =3.1415926535897932384626433832795028841972 ;
static const double p[g+2] = {0.99999999999980993, 676.5203681218851,
-1259.1392167224028, 771.32342877765313, -176.61502916214059,
12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
1.5056327351493116e-7};

dcomplex gamma(dcomplex z)
{
	return z;
}
*/



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

void print_usage()
  /* called when wrong arguments passed via command line */
{
	printf("Usage: double-fourier%s (hur-file) (facs-file) (ofile)\n",VERSION);
	printf("  (hur-file)  - file with hurwitz values\n");
	printf("  (facs-file) - file with factors\n");
	printf("  (ofile)     - output file\n");
//	printf("  (zeta-file) - file with zeta/gamma values\n");
	exit(0);
}

void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
	std::cout << error_string << endl;
  exit(0);
}

// read the factors file
// file starts at 3
// each line consists of primitive root (0 if there isn't one)
// phi(n)
// num_facs = number of factors of the form p^n with p prime
// the factors from lowest p up in the form p then p^n
bool read_factors(factor *factors, unsigned int n, ifstream &facs_file)
{
  unsigned int i,j;
  for(i=3;i<=n;i++)
    {
      facs_file >> factors[i].pr >> factors[i].phi >> factors[i].num_facs;
      for(j=0;j<factors[i].num_facs;j++)
		  facs_file >> factors[i].primes[j] >> factors[i].facs[j];
    }
  return(true);
}

void do_1d_fft(dcomplex *fft_vec, unsigned int num_ffts, unsigned int fft_len,
			   dcomplex *w, dcomplex *a, dcomplex *b, dcomplex *b_spare, dcomplex *b_star_conj)
{
	unsigned int conv_size,i;

    init_bluestein_fft(fft_len,&conv_size,w,b,b_star_conj);

	for(i=0;i<num_ffts;i++)
		bluestein_fft(1,&fft_vec[i*fft_len],fft_len,w,a,b,b_star_conj,conv_size,b_spare); 
}

unsigned int fft_2_count=0;
inline void do_2_fft(unsigned int stride, dcomplex *fft_vec)
{
	dcomplex temp=fft_vec[0];
	fft_vec[0]+=fft_vec[stride];
	fft_vec[stride]=temp-fft_vec[stride];
	fft_2_count++;
}

unsigned int fft_4_count=0;
inline void do_4_fft(unsigned int stride, dcomplex *fft_vec)
{
	dcomplex ac=fft_vec[0]+fft_vec[stride*2];
	dcomplex a_c=fft_vec[0]-fft_vec[stride*2];
	dcomplex bd=fft_vec[stride]+fft_vec[stride*3];
	dcomplex ib_d=fft_vec[stride]-fft_vec[stride*3];
	int_double ib_d_re=ib_d.real;
	fft_vec[0]=ac+bd;
	fft_vec[stride*2]=ac-bd;
	ib_d.real=-ib_d.imag;
	ib_d.imag=ib_d_re;
	fft_vec[stride*3]=a_c+ib_d;
	fft_vec[stride]=a_c-ib_d;
	fft_4_count++;
}

unsigned int fft_power_2_count=0;
inline void power_2_fft(unsigned int n, unsigned int stride, dcomplex *fft_vec)
{
	unsigned int i,j;
	fft_power_2_count++;
	if(stride==1)
		fft(fft_vec,n,w);
	else
	{
		for(i=0,j=0;i<n*stride;i+=stride,j++)
			a[j]=fft_vec[i];
		fft(a,n,w);
		for(i=0,j=0;i<n*stride;i+=stride,j++)
			fft_vec[i]=a[j];
	}
}


// don't call this with n==0!
bool bin_power_2(unsigned int n)
{
	while((n&1)!=1)
		n>>=1;
	return(n==1);
}


void do_nd_ffts(dcomplex *fft_vec, unsigned int num_ffts, unsigned int rank, const int *dims,unsigned int dist,
				dcomplex *w, dcomplex *a, dcomplex *b, dcomplex *b_spare, dcomplex *b_star_conj)
{
	// do <num_ffts> <rank>-d ffts on <fft_vec> each of length <length> where <dist>=prod(dims) apart

	unsigned int i=0,j,k,l,conv_size,stride=dist,depth;

	for(i=0;i<rank;i++)
	{
		if(dims[i]==2)
			{	
				depth=dist/stride;
				stride>>=1;
				for(j=0;j<num_ffts*dist;j+=dist)
					for(k=0;k<depth*(stride<<1);k+=(stride<<1))
						for(l=0;l<stride;l++)
							do_2_fft(stride,&fft_vec[j+k+l]);
				continue;
			}
		if(dims[i]==4)
		{
			depth=dist/stride;
			stride>>=2;
			for(j=0;j<num_ffts*dist;j+=dist)
				for(k=0;k<depth*(stride<<2);k+=(stride<<2))
					for(l=0;l<stride;l++)
						do_4_fft(stride,&fft_vec[j+k+l]);
			continue; 
			}
		
		if((bin_power_2(dims[i])))
		{
			initfft(dims[i],w);
			depth=dist/stride;
			stride/=dims[i];
			for(j=0;j<num_ffts*dist;j+=dist)
				for(k=0;k<depth*dims[i]*stride;k+=dims[i]*stride)
					for(l=0;l<stride;l++)
					    power_2_fft(dims[i],stride,&fft_vec[j+k+l]);
			continue;
		}


			depth=dist/stride;
			stride/=dims[i];
			init_bluestein_fft(dims[i],&conv_size,w,b,b_star_conj);
			for(j=0;j<num_ffts*dist;j+=dist)
				for(k=0;k<depth*dims[i]*stride;k+=dims[i]*stride)
					for(l=0;l<stride;l++)
						bluestein_fft(stride,&fft_vec[j+k+l],dims[i],w,a,b,b_star_conj,conv_size,b_spare);
	}
}

// on exit q_n[(pr^n) mod q]=n
void fill_q_ns (unsigned int *q_n,// unsigned int *q_neg_n,
		unsigned int pr, unsigned int phi_q,
		unsigned int q)
{
  unsigned int i,pr_n;

  q_n[1]=0;
  q_n[pr]=1;
  pr_n=pr;
  for(i=2;i<phi_q;i++)
  {
	  pr_n=(pr_n*pr)%q;
	  q_n[pr_n]=i;
  }
}


void fill_q_ns_2s(unsigned int *q_n, unsigned int q)
{
  unsigned int i,pr;
  pr=1;
  for(i=0;i<(q>>2);i++)
  {
	  q_n[pr]=i;
	  q_n[q-pr]=i+(q>>2);
	  pr=(pr*5)%q;
  }
}

// on exit, offset[i] is where the i'th fraction a/q goes in the fft vector
void make_offsets(unsigned int q, unsigned int *q_n, unsigned int *offset_vec, factor *factors,
				  unsigned int *a_n)
{
	unsigned int i,j,fac,fac1,offset,offset1,ptr;
	
	ptr=0;   // ptr into offset vec
	fac=factors[q].facs[0];
	for(i=0;i<factors[q].phi;i++)
		{
			offset=q_n[a_n[i]%fac];
			offset1=fac;
			for(j=1;j<factors[q].num_facs;j++)
			{
				fac1=factors[q].facs[j];
				offset*=factors[fac1].phi;
				offset+=q_n[a_n[i]%fac1+offset1];
				offset1+=fac1;
			}
			offset_vec[ptr++]=offset;
		}
}

bool non_prim(unsigned int a_n, unsigned int q, factor *factors)
{
	unsigned int i;

	for(i=0;i<factors[q].num_facs;i++)
		if(a_n%factors[q].facs[i]==1)
			return(false);
	return(true);
}

/*
void make_omegas(dcomplex *omegas, dcomplex *fft_vec, unsigned int phi_q, int_double im_s, unsigned int q)
{
	unsigned int i;

	dcomplex q_it_2; // q^(-it/2)
	dcomplex z;

	int_double ln_q;

	for(i=1;i<phi_q;i++) // i=0 is never primitive
	{
		ln_q=log(int_double(q))*(-im_s/2.0);
		sin_cos(ln_q,&q_it_2.imag,&q_it_2.real);
		z=fft_vec[i]*q_it_2;
		omegas[i]=conj(z/sqrt(norm(z)));
	}
}
*/

dcomplex make_omega(dcomplex &z)
{
	dcomplex z1;
	if(contains_zero(z.real)&&contains_zero(z.imag))
	{
		printf("Fatal Error in make_omega. z contains zero. Exiting.\n");
		exit(0);
	}

	z1=conj(z/sqrt(norm(z)));
	// now put z1 in right half plane
	if(z1.real.right>=0.0)
		return(-z1);
	if(z1.real.left>=0.0)
		return(z1);
	// z1 straddles imaginary axis
	printf("Worry in make_omega - z1 straddles imaginary axis.\nz= ");
	print_int_complex(z);printf("\nz1= ");
	print_int_complex(z1);printf("\nUsing z1 value.\n");
	return(z1);
}

inline bool prim_p(unsigned int i, unsigned int p)
{
	return((i%p)!=0);
}

// returns true if the relevant character is odd
// i.e. if chi(-1)=-1
// true if the parity of the co-ordinates is odd
inline bool neg_one_p(unsigned int *coords, unsigned int n)
{
	unsigned int i,sum;
	sum=0;
	for(i=0;i<n;i++)
		sum+=coords[i];
	return(sum&1);
}

void prep_omegas(dcomplex *omegas, unsigned int q, unsigned int phi_q, unsigned int *q_n,
				 dcomplex *w, dcomplex *a, dcomplex *b, dcomplex *b_star_conj,
				 unsigned int conv_size, dcomplex *b_spare)
{
	unsigned int n;
	int_double theta,q_minus_half;

	q_minus_half=exp(log(int_double(q))*(-0.5));
	for(n=1;n<q;n++)
		if(co_prime(n,q))
		{
			theta=d_two_pi*n;
			theta/=q;
			sin_cos(theta,&omegas[q_n[n]].imag,&omegas[q_n[n]].real);
		}
	bluestein_fft(1,omegas,phi_q,w,a,b,b_star_conj,conv_size,b_spare);
	for(n=0;n<phi_q;n++)
		omegas[n]*=q_minus_half;
}

void prep_omegas_nd(unsigned int *offsets,unsigned int q, dcomplex *omegas,unsigned int n_dims,int *dims,
					unsigned int phi_q,dcomplex *w,dcomplex *a,dcomplex *b,dcomplex *b_spare,dcomplex *b_star_conj)
{
	unsigned int n,n1;
	int_double theta,q_minus_half;

	q_minus_half=exp(log(int_double(q))*(-0.5));
	n1=0;
	for(n=1;n<q;n++)
		if(co_prime(n,q))
		{
			theta=d_two_pi*n;
			theta/=q;
			sin_cos(theta,&omegas[offsets[n1]].imag,&omegas[offsets[n1]].real);
			n1++;
		}

	do_nd_ffts(omegas,1,n_dims,dims,phi_q,w,a,b,b_spare,b_star_conj);

	for(n=0;n<phi_q;n++)
		omegas[n]*=q_minus_half;
}



dcomplex finish_omega(dcomplex &omega, bool neg_one) // do the i^(-delta) and sqrt'ing
{
	dcomplex res,one;
	one=c_one;
	if(neg_one) // mult by -i
	{
		res.real=omega.imag;
		res.imag=-omega.real;
	}
	else
		res=omega;

	res=one/sqrt(res);
	if(res.real.right>=0)
	{
		printf("sqrt returned something in left half plane.\n");
		return(-res);
	}
	if(res.real.left>=0)
		return(res);
	printf("Dodgy omega value ");print_int_complex(res);printf("\n");
	return(res);
}

unsigned int num_prims(unsigned int q,factor *factors) // returns number of primitive characters mod q
{
	unsigned int res=1,i,p;
	for(i=0;i<factors[q].num_facs;i++)
	{
		if(factors[q].facs[i]==factors[q].primes[i])
			res*=factors[q].facs[i]-2;
		else
		{
			p=factors[factors[q].facs[i]].phi;
			res*=(p-p/factors[q].primes[i]);
		}
	}
	return(res);
}


void make_l(unsigned int q, 
			unsigned int num_s, 
			factor *factors,
			unsigned int *q_n,
			unsigned int step,
			im_s *im_s_vec,
			unsigned int *a_n,
			unsigned int *offset_vec,
			dcomplex *omegas,
			FILE *hur_file, 
			FILE *out_file)
{
	unsigned int phi_q=factors[q].phi,fac,coords[MAX_FACS],conv_size;
	dcomplex omega,omega_a,z,z1;
	int dims[MAX_FACS];
	int no_dims;
	bool power_2,primitive,neg_one;
	unsigned int pr=factors[q].pr,n_prims;
	unsigned int i,j,fft_ptr,offset,s_done,ffts_done,s_ptr;
	int file_pos;

	printf("Processing Q=%d\n",q);
	file_pos=ftell(hur_file);  // remember where in file we are

	fwrite(&q,sizeof(unsigned int),1,out_file);
	n_prims=num_prims(q,factors); // no of prim characters mod q
	fwrite(&n_prims,sizeof(unsigned int),1,out_file);

	j=0;

//	make_im_s_vec(im_s_vec,q,num_s,gamma_s,gamma_s_a);
	for(i=1;i<q;i++)
		if(co_prime(i,q))
			a_n[j++]=i;

	if(pr)  // q has a primitive root so nice and easy 1-d FFT
	{
		init_bluestein_fft(phi_q,&conv_size,w,b,b_star_conj);
		
		fill_q_ns(q_n,pr,phi_q,q);
		prep_omegas(omegas,q,phi_q,q_n,w,a,b,b_star_conj,conv_size,b_spare);

		fft_ptr=0;
		s_done=0;
		s_ptr=0;
		ffts_done=0;
		while(s_done<num_s)
		{
			while(fft_ptr<FFT_VEC_SIZE-phi_q)   // can fit in another s
			{
				for(i=1;i<q;i++)
					if(co_prime(i,q))
						fread(&fft_vec[fft_ptr+q_n[i]],sizeof(dcomplex),1,hur_file); //Zeta(s,a/q)*q^(-s/2)
				fft_ptr+=phi_q;
				fseek(hur_file,sizeof(dcomplex)*(step-phi_q),SEEK_CUR);  // skip the rest of the fractions
				s_done++;
				ffts_done++;
				if(s_done==num_s)
					break;
			}
			
			for(i=0;i<ffts_done;i++)
				bluestein_fft(1,&fft_vec[i*phi_q],phi_q,w,a,b,b_star_conj,conv_size,b_spare); 

			// now fft_vec contains L(chi,s)*q^(s/2)
			for(i=1;i<phi_q;i++)  // i=0 corresponds to principal chi
				if(prim_p(i,factors[q].primes[0]))				
				{
					//printf("\n");
					fwrite(&i,sizeof(unsigned int),1,out_file);
					fwrite(&ffts_done,sizeof(unsigned int),1,out_file);
					//printf("q: %d i: %d\n",q,i);
					fft_ptr=i;
					neg_one=(i&1);
					
					if(ffts_done==s_done) // first time through, so finish off omega values
						{
							omegas[i]=finish_omega(omegas[i],neg_one);
							fwrite(&omegas[i],sizeof(dcomplex),1,out_file);
							fwrite(&neg_one,sizeof(bool),1,out_file);
						}

					for(j=0;j<ffts_done;j++)
					{
						z=fft_vec[fft_ptr];			

						if(neg_one) // chi(-1)=-1
							z=z*omegas[i]*im_s_vec[j+s_done-ffts_done].lambda_s_a;
						else
							z=z*omegas[i]*im_s_vec[j+s_done-ffts_done].lambda_s;

						
							
						//printf("im_s: %10.8e ",im_s_vec[j].im_s);
						//print_int_complex(z);
						//printf("\n");

						if(!contains_zero(z.imag))
						{
							printf("Imaginary part of z does not contain zero.\n");
						}

						fwrite(&z.real,sizeof(int_double),1,out_file);
						fft_ptr+=phi_q;
					}
				}
			fft_ptr=0;
			ffts_done=0;
		}
	}
	else // q doesn't have a primitive root
	{
		no_dims=factors[q].num_facs;
		fac=factors[q].facs[0];        // the first p^n
		power_2=(factors[fac].pr==0);  // no conductor => p=2, n>=3
		if(power_2)
		{
			no_dims++;
			fill_q_ns_2s(q_n,fac);    // use the {-1,1}X{5} trick
			for(i=1;i<factors[q].num_facs;i++) // move all the dimensions up one
				dims[i+1]=factors[factors[q].facs[i]].phi;
			dims[1]=factors[q].facs[0]/4;
			dims[0]=2;                         // slip in a two
		}
		else
		{
			fill_q_ns(q_n,factors[fac].pr,factors[fac].phi,fac); // use the generator
			for(i=0;i<factors[q].num_facs;i++)
				dims[i]=factors[factors[q].facs[i]].phi;
		}

		offset=fac;
		for(i=1;i<factors[q].num_facs;i++)  // do the rest on the factors, all will have generators
		{
			fac=factors[q].facs[i];			
			pr=factors[fac].pr;
			fill_q_ns(&q_n[offset],pr,factors[fac].phi,fac);  // use the generator
			offset+=fac;
		}
		fft_ptr=0;
		s_done=0;
		ffts_done=0;
		make_offsets(q,q_n,offset_vec,factors,a_n);    // reverse q_n so we know how to populate fft vector
		prep_omegas_nd(offset_vec,q,omegas,no_dims,dims,phi_q,w,a,b,b_spare,b_star_conj);

		

		while(s_done<num_s)
		{
			while(fft_ptr<=FFT_VEC_SIZE-phi_q) // can fit another s
			{
				for(i=0;i<phi_q;i++)
					fread(&fft_vec[fft_ptr+offset_vec[i]],sizeof(dcomplex),1,hur_file);
				fft_ptr+=phi_q;
				fseek(hur_file,sizeof(dcomplex)*(step-phi_q),SEEK_CUR);  // skip the rest of the fractions
				s_done++;
				ffts_done++;
				if(s_done==num_s)
					break;
			} 

			do_nd_ffts(fft_vec,ffts_done,no_dims,dims,phi_q,w,a,b,b_spare,b_star_conj);

			for(i=0;i<factors[q].num_facs;i++)
				coords[i]=0;

			for(i=0;i<phi_q;i++)

			{
				//printf("%10.8e+%10.8ei\n",fft_vec[i].real.left,fft_vec[i].imag.left);
				primitive=true;

				for(j=0;j<factors[q].num_facs;j++)
					if(coords[j]%factors[q].primes[j]==0)
					{
						primitive=false;
						break;
					}
					
				if(primitive)
				{
					fwrite(&i,sizeof(unsigned int),1,out_file);
					fwrite(&ffts_done,sizeof(unsigned int),1,out_file);

					//printf("\n");
					fft_ptr=i;

					
					for(j=0;j<ffts_done;j++)
					{
						z=fft_vec[fft_ptr];

						neg_one=neg_one_p(coords,factors[q].num_facs);
						if(power_2&&(i<(phi_q>>1)))
							neg_one=!neg_one;
							
						if((j==0)&&(ffts_done==s_done)) // first time for this im_s and chi
						{
//							if(i==24411)
//							printf("%d z=[%10.8e,%10.8e]+i[%10.8e,%10.8e]\n",i,z.real.left,-z.real.right,
//								z.imag.left,-z.imag.right);
							omegas[i]=finish_omega(omegas[i],neg_one);
							fwrite(&omegas[i],sizeof(dcomplex),1,out_file);
							fwrite(&neg_one,sizeof(bool),1,out_file);
						}
						
						if(neg_one)
								z*=omegas[i]*im_s_vec[j+s_done-ffts_done].lambda_s_a;
							else
								z*=omegas[i]*im_s_vec[j+s_done-ffts_done].lambda_s;

							
						//printf("%10.8e %20.18e %20.18e\n",im_s_vec[j].im_s,z.real.left,-z.real.right);
						
						if(!contains_zero(z.imag))
						{
							printf("Non Zero Imag Part %d %d %d\n",q,i,j);
							//printf("Imaginary part of z does not contain zero.\n");
							//exit(0);
						}

						fwrite(&z.real,sizeof(int_double),1,out_file);
						fft_ptr+=phi_q;
					}
				}

				j=factors[q].num_facs-1;
				coords[j]++;
				while(coords[j]==factors[factors[q].facs[j]].phi)
				{
					coords[j]=0;
					if(j==0)
						break;
					j--;
					coords[j]++;
				}
			}
			ffts_done=0;
			fft_ptr=0;
		}
	}
	fseek(hur_file,sizeof(dcomplex)*phi_q+file_pos,SEEK_SET);
}

int main(int argc, char **argv)
{
	unsigned int q_start,q_end,step,num_s,q,i,*q_n,*offset_vec,*a_n;
	FILE *hur_file,*out_file;
	ifstream facs_file;
	factor *factors;
	dcomplex *omegas;
	im_s *im_s_vec;
	

	clock_t no_clicks;

	no_clicks=clock(); // start timing

	_fpu_rndd;

	if(argc!=4)
		print_usage();
	hur_file=fopen(argv[1],"rb");
	if(!hur_file)
		fatal_error("Couldn't open hurwitz data file. Exiting.\n");
	facs_file.open(argv[2]);
	if(!facs_file.is_open())
		fatal_error("Couldnt open factors file. Exiting.\n");
	out_file=fopen(argv[3],"wb");
	if(!out_file)
		fatal_error("Couldn't open output file. Exiting.\n");

	fread(&q_start,sizeof(unsigned int),1,hur_file);  // lowest q in file
	fread(&q_end,sizeof(unsigned int),1,hur_file);    // highest q in file
	fread(&step,sizeof(unsigned int),1,hur_file);     // total a/q's in file
	fread(&num_s,sizeof(unsigned int),1,hur_file);    // number of s values in file
	fwrite(&num_s,sizeof(unsigned int),1,out_file);

	printf("Running from q= %d to %d.\n",q_start,q_end);

	if(!(im_s_vec=(im_s *) calloc(num_s,sizeof(im_s))))
		fatal_error("Error allocating memory for im_s_vec. Exiting.\n");
/*	
	if(!(gamma_s=(dcomplex *) calloc(sizeof(dcomplex)*num_s)))
		fatal_error("Error allocating memory for im_s_vec. Exiting.\n");
	
	if(!(gamma_s_a=(dcomplex *) calloc(sizeof(dcomplex)*num_s)))
		fatal_error("Error allocating memory for gamma_s_a. Exiting.\n");
*/

	// read all the imaginary parts of the num_s s values
	// should do this on the fly in make_l
	
	for(i=0;i<num_s;i++)
	{
		fread(&im_s_vec[i].im_s,sizeof(double),1,hur_file);
		fread(&im_s_vec[i].lambda_s,sizeof(dcomplex),1,hur_file);
		fread(&im_s_vec[i].lambda_s_a,sizeof(dcomplex),1,hur_file);
	}
	fwrite(im_s_vec,sizeof(im_s),num_s,out_file);

	if(!(factors=(factor *) calloc(q_end+1,sizeof(factor))))
		fatal_error("Error allocating memory for factors. Exiting.\n");

	if(!read_factors(factors,q_end,facs_file))
		fatal_error("Error reading factor file. Exiting.\n");
	facs_file.close();

	if(!(omegas=(dcomplex *) calloc(q_end-1,sizeof(dcomplex))))
		fatal_error("Error allocating memory for omegas. Exiting.\n");

	if(!(q_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
		fatal_error("Error allocating memory for q_n. Exiting.\n");

	if(!(a_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
		fatal_error("Error allocating memory for a_n. Exiting.\n");

	if(!(offset_vec=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
		fatal_error("Error allocating memory for offset_vec. Exiting.\n");

	for(q=q_start;q<=q_end;q++)
	{
		// hur_file is pointing to H(s_0,1/q_start)
		if((q&3)!=2)
			make_l(q,num_s,factors,q_n,step,
			im_s_vec,offset_vec,a_n,omegas,hur_file,out_file);
	}


	fclose(out_file);
	fclose(hur_file);

	printf("Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
	printf("%ld length 2 FFTs computed.\n",fft_2_count);
	printf("%ld length 4 FFTs computed.\n",fft_4_count);
	printf("%ld power of 2 length FFTs computed.\n",fft_power_2_count);
	return(0);
}



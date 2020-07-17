//
// Version 1.2 computes DFT for 1/2 of Bluestein convolution once only.
//
// Version 1.3 does 1/2 length DFT's and combines.
//

#ifndef INT_FFT
#define INT_FFT

// load Windows memory management stuff
#ifndef LINUX
#include <windows.h>
#endif

//#include "../includes/fft_defs.h"
#define MAX_Q (300000) // actually 750000 for odd q
// the largest phi is 3e5 for even, 1.5e5 for odd
// so the largest DFT is 1.5e5
// and we use the Hermitian trick to halve that
#define MAX_CONV (1<<18)//23) // smallest pwr of 2 >=phi(Q)-1
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


inline long unsigned int gcd (long unsigned int a, long unsigned int b)
/* Euclid algorithm gcd */
{
	long unsigned int c;
	while(a!=0)
	{
		c=a;
		a=b%a;
		b=c;
	};
	return(b);
};

inline int co_prime(long unsigned int a, long unsigned int b)
{
	return(gcd(a,b)==1);
};


// calculate all the sin_cos factors
void init_ws(int_complex *ws)
{
  unsigned int ws_ptr=1; // this is where n=4 starts
  unsigned int n,i;
  int_double two_pi_n;
  
  ws[0]=c_one;
  for(n=4;n<=MAX_CONV;n+=n)
    {
      ws[ws_ptr++]=c_one;
      two_pi_n=int_double(2.0)/n;
      for(i=1;i<n/2;i++)
	{
	  sin_cospi(two_pi_n*i,&ws[ws_ptr].imag,&ws[ws_ptr].real);
	  ws_ptr++;
	}
    }
}

// perform an in place FFT ,n a power of 2
// F[k]=sum_{j=0..N-1} F[j]e(-kj/N)
void fft(int_complex *x,unsigned int n,int_complex *w)
{
  unsigned int i,j,k,l;
  int_complex *p,*xend=x+n,ctemp;

  /* swap each element with one with bit-reversed index */
  for (i=0,l=n>>1;i<l;++i) 
    {
      /* j = bit reversal of i */
      for (k=1,j=0;k<n;k<<=1) 
	{
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
      for (j=0;j<n/2;j+=l,p++) 
	{
	  ctemp=p[k]*w[j];
	  p[k]=p[0]-ctemp;
	  p[0]=p[0]+ctemp;
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


int_complex _ws[MAX_CONV],a[MAX_CONV],bs[MAX_CONV],b_star_conjs[MAX_CONV];
unsigned int conv_sizes[MAX_DIMS],ws_ptr[MAX_DIMS];

int_complex *simple_dft_ws[MAX_SIMPLE_DFT+1];

// do straight forward O(n^2) DFT for small n.
// NB I can predict which ones so could speed this up a bit.
void simple_dft(unsigned int n, unsigned int stride, int_complex *simple_vec)
{
  unsigned int n_2=(n+1)>>1,i,j,ind;
  int_complex w;
  if(!simple_dft_ws[n]) // haven't done one this big before
    {
      //printf("Preparing n=%ld.\n",n);
      if(!(simple_dft_ws[n]=(int_complex *) _aligned_malloc(sizeof(int_complex)*n_2,16)))
	{
	  printf("Couldn't allocate memory for simple_dft_ws[%ld]. Exiting.\n",n);
	  exit(1);
	}
      simple_dft_ws[n][0]=c_one;
      for(i=1;i<n_2;i++)
	{
	  sin_cospi(int_double(2*i)/n,&simple_dft_ws[n][i].imag,&simple_dft_ws[n][i].real);
	  //			print_int_complex_str("w= ",simple_fft_ws[n][i]);
	}
    }

  for(i=0;i<n;i++)
    {
      a[i]=simple_vec[0];
      for(j=1;j<n;j++)
	{
	  ind=(i*j)%n;
	  if(ind>=n_2)
	    {
	      if(n&1)
		w=conj(simple_dft_ws[n][n-ind]);
	      else
		w=-simple_dft_ws[n][ind-n_2];
	    }
	  else
	    w=simple_dft_ws[n][ind];
	  a[i]+=simple_vec[stride*j]*w;
	}
    }

  for(i=0,j=0;i<n;i++,j+=stride)
    {
      simple_vec[j]=a[i];
      //		print_int_complex_str("fft_vec ",simple_vec[i*stride]);
    }
}




/* circular convolution; f and g must be distinct */
void convolve(int_complex *result,int_complex *f,int_complex *g,unsigned int n)
{
	unsigned int i;
	fft(f,n,&_ws[(n>>1)-1]);
	fft(g,n,&_ws[(n>>1)-1]);
	for (i=0;i<n;i++)
		result[i]=f[i]*g[i];
	ifft(result,n,&_ws[(n>>1)-1]);
}

void init_bluestein_fft(unsigned int n, unsigned int *conv_size,
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

	if(n<=MAX_SIMPLE_DFT)
	  return;
	conv_size[0]=1;
	while(conv_size[0]<2*n-1)
		conv_size[0]<<=1;
	if(conv_size[0]>MAX_CONV)
	  {
	    printf("MAX_CONV=%u exceeded in init_bluestein_fft by n=%u. Exiting.\n",MAX_CONV,n);
	    exit(0);
	  }

	pi_over_n=d_one/n;
	b[0]=c_one;
	b_star_conj[0]=c_one;
	for(i=1;i<n;i++)
	{
		ismod=(long long unsigned int)i*(long long unsigned int)i;
		ismod=ismod%(2*n);
		isqrmod=(unsigned int) ismod;

		sin_cospi(pi_over_n*isqrmod,&b[i].imag,&b[i].real);

		b_star_conj[i]=conj(b[i]);
	}
	for(i=1;i<n;i++)
		b[conv_size[0]-i]=b[i];
	for(i=n;i<=(conv_size[0]-n);i++)
		b[i]=c_zero;
	/*
	printf("in init_bluestein_fft.\n");
	int r_err=-1000,a_err=-1000,sum_r_err=0,sum_a_err=0;
	for(long unsigned int j=1;j<n;j++)
	  {
	    int this_r_err=rel_error(b[j].real);
	    int this_a_err=abs_error(b[j].real);
	    sum_r_err+=this_r_err;
	    sum_a_err+=this_a_err;
	    if(this_r_err>r_err) r_err=this_r_err;
	    if(this_a_err>a_err) a_err=this_a_err;
	  }
	printf("worst case relative error 10^{%d} absolute error 10^{%d}\n",r_err,a_err);
	printf("average relative error 10^{%g} absolute error 10^{%g}\n",sum_r_err/(double) n,sum_a_err/(double) n);
	*/
	fft(b,conv_size[0],_ws+(conv_size[0]>>1)-1);
	/*
	printf("in init_bluestein_fft.\n");
	r_err=-1000;a_err=-1000;sum_r_err=0;sum_a_err=0;
	for(long unsigned int j=0;j<conv_size[0];j++)
	  {
	    int this_r_err=rel_error(b[j].real);
	    int this_a_err=abs_error(b[j].real);
	    sum_r_err+=this_r_err;
	    sum_a_err+=this_a_err;
	    if(this_r_err>r_err) r_err=this_r_err;
	    if(this_a_err>a_err) a_err=this_a_err;
	  }
	printf("worst case relative error 10^{%d} absolute error 10^{%d}\n",r_err,a_err);
	printf("average relative error 10^{%g} absolute error 10^{%g}\n",sum_r_err/(double) conv_size[0],sum_a_err/(double) conv_size[0]);
	exit(0);
	*/
}

void summarise_errors(int_complex *data, long unsigned int len, long unsigned int stride, const char *str)
{
  /*
  int rel_err=-999,abs_err=-999,this_rel_err,this_abs_err,sum_rel_err=0,sum_abs_err=0,count=0;
  for(long unsigned int ptr=0,i=0;i<len;i++,ptr+=stride)
    {
      this_rel_err=rel_error(data[ptr].real);
      this_abs_err=abs_error(data[ptr].real);
      if(this_rel_err&this_abs_err)
	{
	  count++;
	  if(this_rel_err>rel_err) rel_err=this_rel_err;
	  if(this_abs_err>abs_err) abs_err=this_abs_err;
	  sum_rel_err+=this_rel_err;
	  sum_abs_err+=this_abs_err;
	}
    }
  printf("Reporting on %s length=%d\n",str,len);
  printf("Worst case real error = %5d (relative) and %5d (absolute).\n",rel_err,abs_err);
  printf("Average    real error = %5g (relative) and %5g (absolute).\n",sum_rel_err/(double) count,sum_abs_err/(double) count);
  rel_err=-999;abs_err=-999;this_rel_err;this_abs_err;sum_rel_err=0;sum_abs_err=0;count=0;
  for(long unsigned int ptr=0,i=0;i<len;i++,ptr+=stride)
    {
      this_rel_err=rel_error(data[ptr].imag);
      this_abs_err=abs_error(data[ptr].imag);
      if(this_rel_err&this_abs_err)
	{
	  count++;
	  if(this_rel_err>rel_err) rel_err=this_rel_err;
	  if(this_abs_err>abs_err) abs_err=this_abs_err;
	  sum_rel_err+=this_rel_err;
	  sum_abs_err+=this_abs_err;
	}
    }
  printf("Worst case imag error = %5d (relative) and %5d (absolute).\n",rel_err,abs_err);
  printf("Average    imag error = %5g (relative) and %5g (absolute).\n",sum_rel_err/(double) count,sum_abs_err/(double) count);
  */
}


void bluestein_fft(unsigned int stride, int_complex *x, unsigned int n,
				   int_complex *a, int_complex *b,
				   int_complex *b_star_conj, unsigned int conv_size)
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
	
	if(n<=MAX_SIMPLE_DFT)
	  {
	    //summarise_errors(x,n,stride,"Pre simple DFT.");
	    simple_dft(n,stride,x);
	    //summarise_errors(x,n,stride,"Post simple DFT.");
	    //printf("\n\n");
	    return;
	  }
	//summarise_errors(x,n,stride,"Pre Bluestein FFT.");
	//printf("Here with conv_size=%d stride=%d n=%d.\n",conv_size,stride,n);
	a[0]=x[0]; 
	for(i=1,j=stride;i<n;i++,j+=stride)
		a[i]=x[j]*b_star_conj[i];
	for(i=n;i<conv_size;i++)
		a[i]=c_zero;
	//for(i=0;i<conv_size;i++) b_spare[i]=b[i];
	//printf("Doing forward FFT.\n");
	//summarise_errors(a,n,1,"Pre Gauss FFT.");
	fft(a,conv_size,_ws+(conv_size>>1)-1);
	//summarise_errors(a,n,1,"Post Gauss FFT.");
	//printf("Forward FFT done.\n");
	for(i=0;i<conv_size;i++)
	  a[i]*=b[i];
	//summarise_errors(a,n,1,"Pre inverse FFT.");
	ifft(a,conv_size,_ws+(conv_size>>1)-1);
	//summarise_errors(a,n,1,"Post inverse FFT.");
	//printf("Backward FFT done.\n");
	//convolve(a,a,b_spare,conv_size); // a<- a*b
	x[0]=a[0];
	for(i=1,j=stride;i<n;i++,j+=stride)
		x[j]=a[i]*b_star_conj[i];
	//summarise_errors(x,n,stride,"Post Bluestein FFT.");printf("\n\n");
}


#define FFT_VEC_SIZE (MAX_Q)  // each one is 32 bytes for double complex 

int_complex *Fft_vec;//[FFT_VEC_SIZE];



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

unsigned int fft_2_count=0;
inline void do_2_fft(unsigned int stride, int_complex *fft_vec)
{
	int_complex temp=fft_vec[0];
	fft_vec[0]+=fft_vec[stride];
	fft_vec[stride]=temp-fft_vec[stride];
	fft_2_count++;
}

unsigned int fft_4_count=0;
inline void do_4_fft(unsigned int stride, int_complex *fft_vec)
{
	int_complex ac=fft_vec[0]+fft_vec[stride*2];
	int_complex a_c=fft_vec[0]-fft_vec[stride*2];
	int_complex bd=fft_vec[stride]+fft_vec[stride*3];
	int_complex ib_d=fft_vec[stride]-fft_vec[stride*3];
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
inline void power_2_fft(unsigned int n, unsigned int stride, int_complex *fft_vec)
{
	unsigned int i,j;
	fft_power_2_count++;
	if(stride==1)
		fft(fft_vec,n,&_ws[(n>>1)-1]);
	else
	{
		for(i=0,j=0;i<n*stride;i+=stride,j++)
			a[j]=fft_vec[i];
		fft(a,n,&_ws[(n>>1)-1]);
		for(i=0,j=0;i<n*stride;i+=stride,j++)
			fft_vec[i]=a[j];
	}
}

// don't call with n==0
inline int bin_power_2(unsigned int n)
{
	while((n&1)!=1)
		n>>=1;
	return(n==1);
}

void do_nd_fft(int_complex *fft_vec, unsigned int rank, const int *dims,unsigned int dist)
{
  // do a <rank>-d fft on <fft_vec> of length <length> where <dist>=prod(dims) apart
  unsigned int i=0,k,l,stride=dist,depth;
  
  for(i=0;i<rank;i++)
    {
      if(dims[i]==2)
	{	
	  depth=dist/stride;
	  stride>>=1;
	  for(k=0;k<depth*(stride<<1);k+=(stride<<1))
	    for(l=0;l<stride;l++)
	      do_2_fft(stride,&fft_vec[k+l]);
	  continue;
	}
      if(dims[i]==4)
	{
	  depth=dist/stride;
	  stride>>=2;
	  for(k=0;k<depth*(stride<<2);k+=(stride<<2))
	    for(l=0;l<stride;l++)
	      do_4_fft(stride,&fft_vec[k+l]);
	  continue;
	}
      
      if((dims[i]>MAX_SIMPLE_DFT)&&(bin_power_2(dims[i])))
	{
	  depth=dist/stride;
	  stride/=dims[i];
	  for(k=0;k<depth*dims[i]*stride;k+=dims[i]*stride)
	    for(l=0;l<stride;l++)
	      power_2_fft(dims[i],stride,&fft_vec[k+l]);
	  continue;
	}
      /* now checked within bluestein_fft      
      if(dims[i]<=MAX_SIMPLE_DFT)
	{
	  depth=dist/stride;
	  stride/=dims[i];
	  for(k=0;k<depth*dims[i]*stride;k+=dims[i]*stride)
	    for(l=0;l<stride;l++)
	      simple_dft(dims[i],stride,&fft_vec[k+l]);
	  continue;
	}
      */
      
      
      depth=dist/stride;
      stride/=dims[i];
      for(k=0;k<depth*dims[i]*stride;k+=dims[i]*stride)
	for(l=0;l<stride;l++)
	  bluestein_fft(stride,&fft_vec[k+l],dims[i],a,&bs[ws_ptr[i]],
			&b_star_conjs[ws_ptr[i]],conv_sizes[i]);
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

void prep_omegas(int_complex *omegas, unsigned int q, unsigned int phi_q, unsigned int *q_n,
				 int_complex *a, int_complex *b, int_complex *b_star_conj,
				 unsigned int conv_size)
{
  unsigned int n;
  int_double theta=int_double(2.0)/q,q_minus_half;
  q_minus_half=d_one/sqrt(int_double(q));//exp(log(int_double(q))*(-0.5));
  for(n=1;n<q;n++)
    if(co_prime(n,q))
      sin_cospi(theta*n,&omegas[q_n[n]].imag,&omegas[q_n[n]].real);
  bluestein_fft(1,omegas,phi_q,a,b,b_star_conj,conv_size);
  for(n=0;n<phi_q;n++)
    omegas[n]*=q_minus_half;
}

void prep_omegas_new(int_complex *omegas, unsigned int q, unsigned int phi_q, unsigned int *q_n,
				 int_complex *a, int_complex *b, int_complex *b_star_conj,
		     unsigned int conv_size, int_complex *twiddles)
{
  unsigned int n,i;
  int_double theta=int_double(2.0)/q,q_minus_half;
  q_minus_half=d_one/sqrt(int_double(q));//exp(log(int_double(q))*(-0.5));
  for(n=1;n<q;n++)
    if(co_prime(n,q))
      sin_cospi(theta*n,&omegas[q_n[n]].imag,&omegas[q_n[n]].real);
  bluestein_fft(2,omegas,phi_q/2,a,b,b_star_conj,conv_size);
  bluestein_fft(2,omegas+1,phi_q/2,a,b,b_star_conj,conv_size);
  //simple_dft(phi_q,1,omegas);
  //for(n=0;n<phi_q;n++)
  //omegas[n]*=q_minus_half;
  for(i=0;i<phi_q/2;i++)
    {
      int_complex omega=omegas[i+i]+twiddles[i]*omegas[i+i+1];
      omegas[i+i+1]=(omegas[i+i]-twiddles[i]*omegas[i+i+1])*q_minus_half;
      omegas[i+i]=omega*q_minus_half;
    }
}

void prep_omegas_nd(unsigned int *offsets,unsigned int q, int_complex *omegas,unsigned int n_dims,int *dims,
					unsigned int phi_q)
{
	unsigned int n,n1,i,w_ptr=0;
	int_double theta=int_double(2.0)/q,q_minus_half;

	q_minus_half=d_one/sqrt(int_double(q));//exp(log(int_double(q))*(-0.5));
	n1=0;
	for(n=1;n<q;n++)
		if(co_prime(n,q))
		{
			sin_cospi(theta*n,&omegas[offsets[n1]].imag,&omegas[offsets[n1]].real);
			n1++;
		}
	for(i=0;i<n_dims;i++)
		{
			ws_ptr[i]=w_ptr;
			if((dims[i]==2)||(dims[i]==4))
				continue;

			if(bin_power_2(dims[i]))
			{
				w_ptr+=dims[i];
				continue;
			}

			if(dims[i]<=MAX_SIMPLE_DFT)
				continue;

			init_bluestein_fft(dims[i],&conv_sizes[i],&bs[w_ptr],&b_star_conjs[w_ptr]);
			w_ptr+=conv_sizes[i];
		}

		do_nd_fft(omegas,n_dims,dims,phi_q);

		for(n=0;n<phi_q;n++)
			omegas[n]*=q_minus_half;
}



int_complex finish_omega(int_complex &omega, bool neg_one) // do the i^(-delta) and sqrt'ing
{
  int_complex res;
  if(neg_one) // mult by -i and conj
    {
      res.real=omega.imag;
      res.imag=omega.real;
    }
  else
    res=conj(omega);
  res=sqrt(res);
  if(res.real.right>=0)
    {
      printf("sqrt returned something in left half plane.\n");
      return(-res);
    }
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

      

#endif

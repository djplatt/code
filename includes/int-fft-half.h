#ifndef INT_FFT
#define INT_FFT

#include "../includes/fft_defs_half.h"

int_complex *_ws,*a,*bs,*b_spare,*b_star_conjs;

// calculate all the sin_cos factors
void init_ws()
{
	unsigned int ws_ptr=1; // this is where n=4 starts
	unsigned int n,i;
	int_double two_pi_n;
	_ws=(int_complex *) malloc(sizeof(int_complex)*MAX_CONV);
	a=(int_complex *) malloc(sizeof(int_complex)*MAX_CONV);
	bs=(int_complex *) malloc(sizeof(int_complex)*MAX_CONV);
	b_spare=(int_complex *) malloc(sizeof(int_complex)*MAX_CONV);
	b_star_conjs=(int_complex *) malloc(sizeof(int_complex)*MAX_CONV);
	if(!(_ws&&a&&bs&&b_spare&&b_star_conjs))
	  {
	    printf("Fatal error allocating memory.\n");
	    exit(0);
	  }
	_ws[0]=c_one;
	for(n=4;n<=MAX_CONV;n+=n)
	{
		_ws[ws_ptr++]=c_one;
		two_pi_n=int_double(2.0)/n;
		for(i=1;i<n/2;i++)
		{
			sin_cospi(two_pi_n*i,&_ws[ws_ptr].imag,&_ws[ws_ptr].real);
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


unsigned int conv_sizes[MAX_DIMS],ws_ptr[MAX_DIMS];

int_complex *simple_dft_ws[MAX_SIMPLE_DFT+1];

// do straight forward O(n^2) DFT for small n.
// NB I can predict which ones so could speed this up a bit.
void simple_dft(unsigned int n, unsigned int stride, int_complex *simple_vec)
{
	unsigned int n_2=n>>1,i,j,ind;
	int_complex w;
	if(!simple_dft_ws[n]) // haven't done one this big before
	{
		//printf("Preparing n=%ld.\n",n);
		if(!(simple_dft_ws[n]=(int_complex *) _aligned_malloc(sizeof(int_complex)*n_2,16)))
		{
			printf("Couldn't allocate memory for simple_dft_ws[%d]. Exiting.\n",n);
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
				w=-simple_dft_ws[n][ind-n_2];
			else
				w=simple_dft_ws[n][ind];
			a[i]+=simple_vec[stride*j]*w;
		}
	}

	for(i=0;i<n;i++)
	{
		simple_vec[i*stride]=a[i];
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

	conv_size[0]=1;
	while(conv_size[0]<2*n-1)
		conv_size[0]<<=1;

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
}

void bluestein_fft(unsigned int stride, int_complex *x, unsigned int n,
				   int_complex *a, int_complex *b,
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

	if(n<=MAX_SIMPLE_DFT)
	  {
	    simple_dft(n,stride,x);
	    return;
	  }

	a[0]=x[0]; 
	for(i=1,j=stride;i<n;i++,j+=stride)
		a[i]=x[j]*b_star_conj[i];
	for(i=n;i<conv_size;i++)
		a[i]=c_zero;
	for(i=0;i<conv_size;i++)
		b_spare[i]=b[i];
	convolve(a,a,b_spare,conv_size); // a<- a*b
	x[0]=a[0];
	for(i=1,j=stride;i<n;i++,j+=stride)
		x[j]=a[i]*b_star_conj[i];
}


#define FFT_VEC_SIZE (MAX_Q)  // each one is 32 bytes for double complex 

int_complex Fft_vec[FFT_VEC_SIZE],Fft_copy[FFT_VEC_SIZE],Fft_copy1[FFT_VEC_SIZE];



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
  /*
  printf("in do_nd_fft with rank=%ld, dist=%ld and dims=%ld",rank,dist,dims[0]);
  for(int i=1;i<rank;i++)
    printf(",%ld",dims[i]);
  printf("\n");
  */
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
			&b_star_conjs[ws_ptr[i]],conv_sizes[i],b_spare);
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
				 unsigned int conv_size, int_complex *b_spare)
{
  unsigned int n;
  int_double theta=int_double(2.0)/q,q_minus_half;
  q_minus_half=d_one/sqrt(int_double(q));//exp(log(int_double(q))*(-0.5));
  for(n=1;n<q;n++)
    if(co_prime(n,q))
      sin_cospi(theta*n,&omegas[q_n[n]].imag,&omegas[q_n[n]].real);
  bluestein_fft(1,omegas,phi_q,a,b,b_star_conj,conv_size,b_spare);
  for(n=0;n<phi_q;n++)
    omegas[n]*=q_minus_half;
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
      //printf("initialising bluestein dim %ld\n",i);
      init_bluestein_fft(dims[i],&conv_sizes[i],&bs[w_ptr],&b_star_conjs[w_ptr]);
      //printf("bluestein initialised.\n");
      w_ptr+=conv_sizes[i];
    }
  //printf("doing %ld dimension fft\n",n_dims);
  do_nd_fft(omegas,n_dims,dims,phi_q);
  //printf("%ld dimension fft finished\n",n_dims);  

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

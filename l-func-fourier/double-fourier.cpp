/*

File: double-fourier.cpp

Created: 15th September 2008

Version: <v> = 1.0

Last Modified: 15th September 2008

Dialect: C++

Requires: fttw3

Implementation notes: 

Build instructions: g++ -odouble-fourier double-fourier.cpp -lfftw3 -O3

Reads Double Precision Hurwitz Values and performs FFT to calculate
corresponding l-func values

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define _CRT_SECURE_NO_WARNINGS  // shut Microsoft up
#define VERSION "1.0"
#define MAX_FACS (6)  // 2*4*3*5*7*11*13 > 100,000
#define FFT_VEC_SIZE (1<<25)  // each one is 16 bytes for double complex 

#define for_real
#define TWO_GB // tell code to honour 2GB file limit

#ifdef TWO_GB
#define _int64 fpos_t
#define _fseeki64(x,y,z) fseek(x,y,z)
#define _ftelli64(x) ftell(x)
#endif




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <fftw3.h>
#include <time.h>

#define debug printf("got to line %d\n",__LINE__);

using namespace std;

typedef complex<double> dcomplex;

// structure to hold factor information
typedef struct
{
  unsigned int pr;   // = 0 if no primitive root
  unsigned int phi;  // Euler phi(q)
  unsigned int num_facs;  // number of prime power factors
  unsigned int facs[MAX_FACS];  // the factors p^n
  unsigned int primes[MAX_FACS]; // the prime factors p
} factor;

// structure to hold gamma values etc.
typedef struct
{
	double im_s; // imaginary part of s
	dcomplex lambda_s; // q^(s/2)*pi^(-s/2)*gamma(s/2)
	dcomplex lambda_s_a; // q^(s/2)*pi^(-(s+1)/2)*gamma((s+1)/2)
} im_s;

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

	if ( real(z)<0.5 ) 
	{
		return pi / (sin(pi*z)*gamma(1.0-z));
	}
	z -= 1.0;
	dcomplex x=p[0];
	for (int i=1; i<g+2; i++) 
	{
		x += p[i]/(z+dcomplex(i,0));
	}
	dcomplex t = z + (g + 0.5);
	return sqrt(2*pi) * pow(t,z+0.5) * exp(-t) * x;
}

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


// return x^z
dcomplex pow (double x, dcomplex z)
{
	double tlnx,xtosigma;

	tlnx=log(x)*imag(z);
	xtosigma=pow(x,real(z));
	return(dcomplex(xtosigma*cos(tlnx),xtosigma*sin(tlnx)));
}


void print_usage()
  /* called when wrong arguments passed via command line */
{
	printf("Usage: double-fourier%s (hur-file) (facs-file) (ofile)\n",VERSION);
	printf("  (hur-file)  - file with hurwitz values\n");
	printf("  (facs-file) - file with factors\n");
	printf("  (ofile)     - output file\n");
	exit(0);
}

void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
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

static const double two_pi=2*pi;


void do_1d_fft(dcomplex *fft_vec, unsigned int num_ffts, const int fft_len)
{
	fftw_plan p;
#ifdef for_real
	p=fftw_plan_many_dft
	  (1,&fft_len,num_ffts,
	   reinterpret_cast<fftw_complex*>(fft_vec),NULL,1,fft_len,
	   reinterpret_cast<fftw_complex*>(fft_vec),NULL,1,fft_len,
	   FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
#endif
}

void do_nd_ffts(dcomplex *fft_vec, unsigned int num_ffts, unsigned int rank, const int *dims,unsigned int dist)
{
	// do <num_ffts> ffts on <fft_vec> each of length <length> where elements are <dist> apart
	fftw_plan p;
#ifdef for_real
	p=fftw_plan_many_dft(rank,dims,num_ffts,
		reinterpret_cast<fftw_complex*>(fft_vec),NULL,1,dist,
		reinterpret_cast<fftw_complex*>(fft_vec),NULL,1,dist,FFTW_BACKWARD,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
#endif
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

void make_im_s_vec(im_s *im_s_vec,unsigned int q, unsigned int num_s)
{
	unsigned int i;
	dcomplex s,s_2,s_a_2,q_s;
	for(i=0;i<num_s;i++)
	{
		s=dcomplex(0.5,im_s_vec[i].im_s);
		s_2=s/dcomplex(2,0);
		s_a_2=s_2+dcomplex(0.5,0);
		q_s=pow(q,s_2);
		im_s_vec[i].lambda_s=q_s*pow(pi,-s_2)*gamma(s_2);
		im_s_vec[i].lambda_s/=sqrt(norm(im_s_vec[i].lambda_s));
		im_s_vec[i].lambda_s_a=q_s*pow(pi,-s_a_2)*gamma(s_a_2);
		im_s_vec[i].lambda_s_a/=sqrt(norm(im_s_vec[i].lambda_s_a));
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



void make_omegas(dcomplex *omegas, dcomplex *fft_vec, unsigned int phi_q, im_s *im_s_vec)
{
	unsigned int i;
	dcomplex z;

	for(i=0;i<phi_q;i++)
		omegas[i]=sqrt(conj(fft_vec[i])/fft_vec[i]);
}

inline bool prim_p(unsigned int i, unsigned int p)
{
	return((i%p)!=0);
}

inline bool neg_one_p(unsigned int *coords, unsigned int n)
{
	unsigned int i,sum;
	sum=0;
	for(i=0;i<n;i++)
		sum+=coords[i];
	return(sum&1);
}


void make_l(unsigned int q, 
			unsigned int num_s, 
			factor *factors,
			unsigned int *q_n,
			dcomplex *fft_vec,
			unsigned int step,
			im_s *im_s_vec,
			unsigned int *a_n,
			unsigned int *offset_vec,
			dcomplex *omegas,
			FILE *hur_file, 
			FILE *out_file)
{
	unsigned int phi_q=factors[q].phi,fac,coords[MAX_FACS];
	dcomplex omega,omega_a;
	int dims[MAX_FACS];
	bool power_2,primitive,neg_one;
	unsigned int pr=factors[q].pr;
	unsigned int i,j,fft_ptr,offset,s_done,ffts_done,s_ptr;
	dcomplex z;
	double x;
	_int64 file_pos;

	file_pos=_ftelli64(hur_file);  // remember where in file we are

	j=0;

	make_im_s_vec(im_s_vec,q,num_s);
	for(i=1;i<q;i++)
		if(co_prime(i,q))
			a_n[j++]=i;

	if(pr)  // q has a primitive root so nice and easy 1-d FFT
	{
		fill_q_ns(q_n,pr,phi_q,q);
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
						fread(&fft_vec[fft_ptr+q_n[i]],sizeof(dcomplex),1,hur_file);
				fft_ptr+=phi_q;
				_fseeki64(hur_file,
					sizeof(dcomplex)*(step-phi_q)+sizeof(double),
					SEEK_CUR);  // skip the rest of the fractions
				// plus the im s value
				s_done++;
				ffts_done++;
				if(s_done==num_s)
					break;
			}

			do_1d_fft(fft_vec,ffts_done,phi_q);
			if(ffts_done==s_done) // first time through for this q
				make_omegas(omegas,fft_vec,phi_q,im_s_vec);
			for(i=1;i<phi_q;i++)  // i=0 corresponds to principal chi
				if(prim_p(i,factors[q].primes[0]))				
				{
					//printf("\n");
					fwrite(&q,sizeof(unsigned int),1,out_file);
					fwrite(&i,sizeof(unsigned int),1,out_file);
					fft_ptr=i;
					for(j=0;j<ffts_done;j++)
					{
						z=fft_vec[fft_ptr];
						if(i&1) // chi(-1)=-1
							z=z*omegas[i]*im_s_vec[j+s_done-ffts_done].lambda_s_a;
						else
							z=z*omegas[i]*im_s_vec[j+s_done-ffts_done].lambda_s;
						x=real(z);
						fwrite(&x,sizeof(double),1,out_file);
						//printf("j:%d q:%d %d %20.18e+%20.18ei\n",j,q,i,real(z),imag(z));
						fft_ptr+=phi_q;
					}
				}
			fft_ptr=0;
			ffts_done=0;
		}
//		fsetpos(hur_file,&file_pos);   // back at H(s_n,1/q_start)
//		fseek(hur_file,sizeof(dcomplex)*phi_q+sizeof(double),SEEK_CUR); // at H(s_{n+1},1/q_start) 
	}
	else // q doesn't have a primitive root
	{
		fac=factors[q].facs[0];        // the first p^n
		power_2=(factors[fac].pr==0);  // no conductor => p=2, n>=3
		if(power_2)
			fill_q_ns_2s(q_n,fac);    // use the {-1,1}X{5} trick
		else
			fill_q_ns(q_n,factors[fac].pr,factors[fac].phi,fac); // use the generator

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

		while(s_done<num_s)
		{
			while(fft_ptr<FFT_VEC_SIZE-phi_q) // can fit another s
			{
				for(i=0;i<phi_q;i++)
					fread(&fft_vec[fft_ptr+offset_vec[i]],sizeof(dcomplex),1,hur_file);
				fft_ptr+=phi_q;
				_fseeki64(hur_file,
					sizeof(dcomplex)*(step-phi_q)+sizeof(double),
					SEEK_CUR);  // skip the rest of the fractions
				// plus the im s value
				s_done++;
				ffts_done++;
				if(s_done==num_s)
					break;
			}

			if(power_2) // first factor is 2^n, n>2
			{
				for(i=1;i<factors[q].num_facs;i++)
					dims[i+1]=factors[factors[q].facs[i]].phi;
				dims[1]=factors[q].facs[0]/4;
				dims[0]=2;
				do_nd_ffts(fft_vec,ffts_done,factors[q].num_facs+1,dims,phi_q);
			}
			else
			{
				for(i=0;i<factors[q].num_facs;i++)
					dims[i]=factors[factors[q].facs[i]].phi;
				do_nd_ffts(fft_vec,ffts_done,factors[q].num_facs,dims,phi_q);
			}
			if(ffts_done==s_done) // first time through for this q
				make_omegas(omegas,fft_vec,phi_q,im_s_vec); // calculate omega for each character mod q
			for(i=0;i<factors[q].num_facs;i++)
				coords[i]=0;
			for(i=0;i<phi_q;i++)
			{
				primitive=true;
				for(j=0;j<factors[q].num_facs;j++)
					if(coords[j]%factors[q].primes[j]==0)
					{
						primitive=false;
						break;
					}
				if(primitive)
				{
					fwrite(&q,sizeof(unsigned int),1,out_file);
					fwrite(&i,sizeof(unsigned int),1,out_file);

					//printf("\n");
					fft_ptr=i;
					for(j=0;j<ffts_done;j++)
					{
						z=fft_vec[fft_ptr];
						neg_one=neg_one_p(coords,factors[q].num_facs);
						if(power_2&&(i<q>>2))
							neg_one=!neg_one;
						if(neg_one)
							z=z*omegas[i]*im_s_vec[j+s_done-ffts_done].lambda_s_a;
						else
							z=z*omegas[i]*im_s_vec[j+s_done-ffts_done].lambda_s;
						x=real(z);
						fwrite(&x,sizeof(double),1,out_file);
						//printf("aq:%d %d %20.18e+%20.18ei\n",q,i,real(z),imag(z));
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
	_fseeki64(hur_file,sizeof(dcomplex)*phi_q+file_pos,SEEK_SET);
}

int main(int argc, char **argv)
{
	unsigned int q_start,q_end,step,num_s,q,i,*q_n,*offset_vec,*a_n;
	FILE *hur_file,*out_file;
	ifstream facs_file;
	_int64 file_pos;
	factor *factors;
	dcomplex *fft_vec;
	dcomplex *omegas;
	im_s *im_s_vec;
	clock_t no_clicks;

	no_clicks=clock(); // start timing

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

	if(!(im_s_vec=(im_s *) malloc(sizeof(im_s)*num_s)))
		fatal_error("Error allocating memory for im_s_vec. Exiting.\n");

	// read all the imaginary parts of the num_s s values
	// should do this on the fly in make_l
	fread(&im_s_vec[0],sizeof(double),1,hur_file);
	fwrite(&im_s_vec[0],sizeof(double),1,out_file);
	file_pos=_ftelli64(hur_file);
	for(i=1;i<num_s;i++)
	{
		_fseeki64(hur_file,sizeof(dcomplex)*step,SEEK_CUR);
		fread(&im_s_vec[i],sizeof(double),1,hur_file);
		fwrite(&im_s_vec[i],sizeof(double),1,out_file);
	}

	// go back to first a/q
	_fseeki64(hur_file,file_pos,SEEK_SET);

	if(!(factors=(factor *) malloc(sizeof(factor)*(q_end+1))))
		fatal_error("Error allocating memory for factors. Exiting.\n");

	if(!read_factors(factors,q_end,facs_file))
		fatal_error("Error reading factor file. Exiting.\n");
	facs_file.close();

	if(!(fft_vec=(dcomplex *) malloc(sizeof(dcomplex)*FFT_VEC_SIZE)))
		fatal_error("Error allocating memory for fft_vec. Exiting.\n");

	if(!(omegas=(dcomplex *) malloc(sizeof(dcomplex)*(q_end-1))))
		fatal_error("Error allocating memory for omegas. Exiting.\n");

	if(!(q_n=(unsigned int *) malloc(sizeof(unsigned int)*(q_end-1))))
		fatal_error("Error allocating memory for q_n. Exiting.\n");

	if(!(a_n=(unsigned int *) malloc(sizeof(unsigned int)*(q_end-1))))
		fatal_error("Error allocating memory for a_n. Exiting.\n");


	if(!(offset_vec=(unsigned int *) malloc(sizeof(unsigned int)*(q_end-1))))
		fatal_error("Error allocating memory for offset_vec. Exiting.\n");


	for(q=q_start;q<=q_end;q++)
	{
		// hur_file is pointing to H(0.5,1/q_start)
		if((q&3)!=2)
			make_l(q,num_s,factors,q_n,fft_vec,step,
			im_s_vec,offset_vec,a_n,omegas,hur_file,out_file);
	};


	fclose(out_file);
	fclose(hur_file);

	printf("Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
	return(0);
};



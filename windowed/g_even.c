//
// g_even.c
//
// Windowed FFT based L-function calculator
//
//
// Vesrion 1.0 Initial implementation
//
// Last Modified: 14 july 2010
//
// DJ Platt
// University of Bristol
//

#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../includes/mpfi_c.h"
#include "../includes/mpfi_fft.h"
#include "../includes/fft_defs.h"


#define TRUE (0==0)
#define FALSE (1==0)
#define FAILURE (1)
#define debug printf("Reached line number %d\n",__LINE__)
#define bool int
#define MAX_N (24) // largest length of FFT array = 2^(MAX_N)

#define one_over_A ((double) 5.0/64.0)

void print_usage()
{
  printf("Usage:- f_hat_even prec q N h t0 facs_file outfile.\n");
  exit(1);
}

int read_factors(factor *factors, unsigned int n, FILE *facs_file)
{
  unsigned int i,j;
  for(i=3;i<=n;i++)
    {
      fscanf(facs_file,"%d %d %d",&factors[i].pr,&factors[i].phi,&factors[i].num_facs);
      for(j=0;j<factors[i].num_facs;j++)
	fscanf(facs_file,"%d %d",&factors[i].primes[j],&factors[i].facs[j]);
    }
  return(TRUE);
}

inline int ln2(int q)
{
  int res=0;
  while(q)
    {
      q>>=1;
      res++;
    }
  return(res);
}

inline int power_2_p(unsigned int i)
{
  while((i&1)==0)
    i>>=1;
  return(i==1);
}

int *co_primes;

void make_chis(unsigned int q, unsigned int index, factor * factors, mpfi_c_t *chis)
{
  unsigned int fac_ptr,fac,phi,pr,prn,i,pw;
  mpfi_t theta,theta1;
  mpfi_c_t chi,chi_neg;

  mpfi_init(theta);
  mpfi_init(theta1);
  mpfi_c_init(chi);
  mpfi_c_init(chi_neg);

  for(i=1;i<q;i++)
    if(co_primes[i])
      mpfi_c_set(chis[i],c_one);

  for(fac_ptr=0;fac_ptr<factors[q].num_facs;fac_ptr++)
    {
      fac=factors[q].facs[fac_ptr];
      phi=factors[fac].phi;
      pr=factors[fac].pr;

      if(pr==0)
      {
	  mpfi_div_ui(theta,mpfi_2_pi,phi/2);
	  mpfi_mul_ui(theta,theta,index);
	  mpfi_neg(theta,theta);
	  //mpfi_print(theta);
	  prn=5;
	  if((index/(phi/2))&1)
	    //if((index<factors[q].phi/4)||((index>=factors[q].phi/2)&&(index<3*factors[q].phi/4)))
	    for(i=fac-1;i<q;i+=fac)
	      mpfi_c_neg(chis[i],chis[i]);
	  for(pw=1;pw<phi/2;pw++)
	    {
	      mpfi_mul_ui(theta1,theta,pw);
	      mpfi_sin(chi->im,theta1);
	      mpfi_cos(chi->re,theta1);
	      for(i=prn;i<q;i+=fac)
		if(co_primes[i])
		  mpfi_c_mul(chis[i],chis[i],chi);
	      for(i=fac-prn;i<q;i+=fac)
	      {
		  if(co_primes[i])
		      mpfi_c_mul(chis[i],chis[i],chi);
		  if((index/(phi/2))&1)
		  //if((index<factors[q].phi/4)||((index>=factors[q].phi/2)&&(index<3*factors[q].phi/4)))
		    //if(index>=factors[q].phi/2)
		      mpfi_c_neg(chis[i],chis[i]);
	      }
	      prn=(prn*5)%fac;
	    }
	}
      else
	{
	  mpfi_div_ui(theta,mpfi_2_pi,phi);
	  mpfi_mul_ui(theta,theta,index);
	  mpfi_neg(theta,theta);
	  for(prn=pr,pw=1;pw<phi;prn=(prn*pr)%fac,pw++)
	    {
	      mpfi_mul_ui(theta1,theta,pw);
	      mpfi_sin(chi->im,theta1);
	      mpfi_cos(chi->re,theta1);
	      //printf("chi=");mpfi_c_print(chi);
	      for(i=prn;i<q;i+=fac)
		if(co_primes[i])
		  {
		    mpfi_c_mul(chis[i],chis[i],chi);
		  }
	    }
	}
      index/=phi;
    }
  mpfi_clear(theta);
  mpfi_clear(theta1);
  mpfi_c_clear(chi);
  mpfi_c_clear(chi_neg);
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

unsigned int conj_j (unsigned int j, unsigned int q, factor *factors)
{
    unsigned int num_chi=num_prims(q,factors);
    if(q&7) // q not divisible by 8
	return(num_chi-j-1);
    if(j<(num_chi>>1))
	return((num_chi>>1)-j-1);
    return(num_chi+(num_chi>>1)-j-1);
}

bool primitive(unsigned int q, unsigned int index, factor *factors)
{
  unsigned int f;
  if((factors[q].primes[0]==2)&&(factors[q].facs[0]>4))
  {
      if(!(index&1))
	  return(FALSE);
      //index/=2;
      //if((index%(factors[q].facs[0]/4))==0)
      //return(FALSE);
      //index/=(factors[q].facs[0]/4);
  }
  else
  {
      if((index%(factors[factors[q].facs[0]].phi))==0)
	  return(FALSE);
  }
  index/=factors[factors[q].facs[0]].phi;

  for(f=1;f<factors[q].num_facs;f++)
    {
      if(((index%factors[factors[q].facs[f]].phi)%factors[q].primes[f])==0)
	return(FALSE);
      index/=factors[factors[q].facs[f]].phi;
    }

  return(TRUE);
}

void calc_exps(unsigned int q, mpfi_c_t *exps)
{
    unsigned int n;
    mpfi_t theta,n_theta;
    mpfi_init(theta);mpfi_init(n_theta);
    mpfi_div_ui(theta,mpfi_2_pi,q);
    for(n=1;n<q;n++)
    {
	if(!co_primes[n])
	    continue;
	mpfi_mul_ui(n_theta,theta,n);
	mpfi_sin(exps[n]->im,n_theta);
	mpfi_cos(exps[n]->re,n_theta);
    }
    mpfi_clear(n_theta);
    mpfi_clear(theta);
}

mpfi_c_t omega_term;

void calc_omega(mpfi_c_ptr omega, unsigned int q, mpfi_c_t *chis, mpfi_c_t *exps, mpfi_ptr sqrt_q)
{
    unsigned int n;
    mpfi_c_zero(omega);
    for(n=1;n<q;n++)
	if(co_primes[n])
	{
	    //printf("chi[n]=");mpfi_c_print(chis[n]);
	    mpfi_c_mul(omega_term,chis[n],exps[n]);
	    //printf("this term (n=%d)=",n);mpfi_c_print(omega_term);
	    mpfi_c_add(omega,omega,omega_term);
	    //printf("running total=");mpfi_c_print(omega);
	}
    mpfi_c_div_i(omega,omega,sqrt_q);
    //printf("omega before sqrt and /i=");mpfi_c_print(omega);
    mpfi_c_add_ui(omega_term,chis[q-1],1);
    if(mpfi_c_contains_zero(omega_term)) // chi[q-1]=-1
    {
	mpfi_neg(omega->re,omega->re);
	mpfi_swap(omega->re,omega->im);
    }
    mpfi_c_sqrt(omega,omega);
}

bool even_p(unsigned int q, mpfi_c_t *chis)
{
   return((mpfi_is_inside_d(-1.0,chis[q-1]->re))&&(mpfi_has_zero(chis[q-1]->im)));
}

void g_e_twiddle_err(mpfi_c_ptr res, double h, double B)
{
  mpfi_t temp;
  //mpfi_c_set_ui(res,0,0);
  //return;
  mpfi_init(temp);
  //printf("Approximating g_e_twiddle error.\n");
  // need to add 4h/B*sqrt(pi/2)*erfc(B/2/sqrt(2)/h)
  mpfi_set_d(res->im,B);
  mpfi_div_d(res->im,res->im,h); // B/h
  mpfi_sqr(res->im,res->im);
  mpfi_neg(res->im,res->im); // -B^2/h^2
  mpfi_div_ui(res->re,res->im,8); // -B^2/8/h^2
  mpfi_exp(res->re,res->re);
  mpfi_mul_ui(res->re,res->re,8); // 8exp(-B^2/8/h^2)
  printf("1st term=");mpfi_printn(res->re,10);

  mpfi_set_ui(res->im,2);
  mpfi_sqrt(res->im,res->im);
  mpfi_set_d(temp,B/2.0);
  mpfi_div(temp,temp,res->im);
  mpfi_div_d(temp,temp,h);
  mpfi_erfc(temp,temp);
  printf("erfc=");mpfi_printn(temp,10);

  mpfi_div_d(temp,temp,B/2.0);
  mpfi_mul_d(temp,temp,h);
  mpfi_set(res->im,mpfi_2_pi);
  mpfi_sqrt(res->im,res->im);
  mpfi_mul(temp,temp,res->im);
  printf("erfc term=");mpfi_printn(temp,10);
  printf("res=");mpfi_c_printn(res,10);
  mpfi_inc(res->re,temp);
  mpfi_neg(res->im,res->re);
  mpfi_put(res->re,res->im);
  mpfi_set(res->im,res->re);
  mpfi_clear(temp);
}

void g_e1(mpfi_c_ptr res, mpfi_c_ptr temp, double t, double t0, double h, unsigned int k)
{
  double t1=t+t0;
  //printf("temp is ");mpfi_c_printn(temp,10);
  mpfi_set_d(temp->re,0.25);
  mpfi_set_d(temp->im,t1/2.0);
  //printf("Taking lng of ");mpfi_c_printn(temp,10);
  mpfi_c_lngc(res,temp);
  mpfi_set_d(temp->im,t1);
  mpfi_mul(temp->im,temp->im,mpfi_pi_by_4);
  mpfi_add(res->re,res->re,temp->im);
  mpfi_set_d(temp->im,t);
  mpfi_div_d(temp->im,temp->im,h);
  mpfi_sqr(temp->im,temp->im);
  mpfi_div_ui(temp->im,temp->im,2);
  mpfi_sub(res->re,res->re,temp->im);
  if((t!=0.0)&&(k!=0))
    {
      mpfi_set_ui(temp->re,0);
      mpfi_mul_d(temp->im,mpfi_2_pi,t);
      mpfi_neg(temp->im,temp->im);
      mpfi_c_log(temp,temp);
      mpfi_c_mul_ui(temp,temp,k); //k log(-2 pi t)
      //printf("Adding ");mpfi_c_printn(temp,10);
      mpfi_c_add(res,res,temp);
    }
  mpfi_c_exp(res,res);
}

  
//
// set Fft_vec=g_e_twiddle(m/A;k) for m=-N/2..N/2-1
// g_e_twiddle(m/A)=sum g_e(m/A+lB;k) for l in Z
//
void g_e(double t0, double h, double B, unsigned int N, unsigned int k, mpfi_c_t *Fft_vec)
{
  int i;
  double t=-one_over_A*(N>>1);
  mpfi_c_t z;
  //mpfi_c_init(z);
  mpfi_c_init(z);
  //mpfi_set_d(z->re,0.25);
  for(i=0;i<N;i++,t+=one_over_A)
    g_e1(Fft_vec[i],z,t,t0,h,k);
  mpfi_c_clear(z);
}

void g_e_twiddle(double h, double B, unsigned int N, mpfi_c_t *Fft_vec)
{
  unsigned int i;
  mpfi_c_t(err);
  mpfi_c_init(err);
  g_e_twiddle_err(err,h,B);
  printf("g_e_twiddle error=");mpfi_c_printn(err,10);
  for(i=0;i<N;i++)
    mpfi_c_inc(Fft_vec[i],err);
  mpfi_c_clear(err);
}

void G_e_twiddle_err(mpfi_c_ptr err,mpfi_ptr A, double h, double t0)
{
  mpfi_t x;
  mpfi_init(x);
  mpfi_mul(err->im,A,mpfi_pi);
  mpfi_set_d(x,1.0);
  mpfi_div(x,x,err->im);
  mpfi_add_ui(x,x,1); // 1+1/Pi/A
  mpfi_sub_d(err->im,err->im,0.5); //Pi*a-0.5
  mpfi_set_d(err->re,1.0);
  mpfi_div_d(err->re,err->re,h);
  mpfi_sqr(err->re,err->re);
  mpfi_mul_d(err->re,err->re,0.125); // 1/8/h^2
  mpfi_sub(err->im,err->re,err->im); // 1/8/h^2+1/2-Pi*A
  mpfi_exp(err->im,err->im);
  mpfi_mul(err->im,err->im,x);
  mpfi_mul(err->im,err->im,mpfi_2_pi);
  mpfi_mul_d(err->im,err->im,h);
  mpfi_neg(x,err->im);
  mpfi_put(err->im,x);
  mpfi_set(err->re,err->im);
  mpfi_clear(x);
}

//
// calculate G_e_twiddle = -1/A * iDFT (g_e_twiddle)
//
void G_e_twiddle(unsigned int N, mpfi_c_t *ws, mpfi_c_t *Fft_vec)
{
  int i;
  nifft(Fft_vec,N,ws); // do iFFT without division by N
  for(i=0;i<N;i++)
    if(i&1)
      mpfi_c_mul_d(Fft_vec[i],Fft_vec[i],one_over_A);
    else
      mpfi_c_mul_d(Fft_vec[i],Fft_vec[i],-one_over_A);

}

//
// calculate G_e = G_e_twiddle+error
//
void G_e(mpfi_ptr A, double h, double t0, unsigned int N, mpfi_c_t *Fft_vec)
{ 
  mpfi_c_t err;
  unsigned int i;
  mpfi_c_init(err);
  G_e_twiddle_err(err,A,h,t0);
  printf("G_e_twiddle_err=");mpfi_c_printn(err,10);
  for(i=0;i<N;i++) 
      mpfi_c_inc(Fft_vec[i],err);
  mpfi_c_clear(err);
}

int main(int argc, char **argv)
{
  int q1,prec,N1;
  unsigned int q,i,N;
  double h,t0,B;
  FILE *outfile,*facs_file;
  factor *factors;
  mpfi_c_t *chis,*exps,omega,*Fft_vec,*ws;
  mpfi_t sqrt_q,A;

  if(argc!=8)
    print_usage();

  prec=atoi(argv[1]);
  if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    print_usage();
  printf("Running at %d bits of precision.\n",prec);
    
  mpfi_c_setup(prec);
  mpfi_c_set_ui(ln_gamma_err1,0,0);printf("Log Gamma Error set to ");mpfi_c_printn(ln_gamma_err1,10);

  mpfi_c_init(omega);
  mpfi_c_init(omega_term);

  q1=atoi(argv[2]);
  if((q1<3)||((q1&3)==2))
    print_usage();
  //printf("q=%d\n",q1);
  q=(unsigned int) q1;

  if(!(co_primes=(int *) malloc (sizeof(int)*q)))
    {
      printf("Failed to allocate memory for co_primes. Exiting\n");
      exit(1);
    }
  co_primes[0]=FALSE;
  co_primes[1]=TRUE;
  for(i=2;i<q;i++)
    co_primes[i]=co_prime(i,q);

  N1=atoi(argv[3]);
  if((N1<=0)||(N1>MAX_N))
    print_usage();
  N=1;	
  for(i=0;i<N1;i++)
    N<<=1;
  printf("N set to %d\n",N);
  mpfi_init(A);
  mpfi_set_ui(A,1);
  mpfi_div_d(A,A,one_over_A);
  printf("A set to ");mpfi_printn(A,10);
  B=N*one_over_A;
  printf("B set to %f\n",B);
  
  h=atof(argv[4]);
  if(h<=0.0)
    print_usage();
  printf("h set to %f\n",h);

  t0=atof(argv[5]);
  if(t0<=0.0)
    print_usage();
  printf("t0 set to %f\n",t0);
  
  facs_file=fopen(argv[6],"r");
  if(!facs_file)
    {
      printf("Failed to open file %s for input. Exiting.\n",argv[3]);
      exit(FAILURE);
    }
  
  factors=(factor *) malloc(sizeof(factor)*(q+1));
  if(!factors)
    {
      printf("Fatal error allocating memory for factors. Exiting.\n");
      exit(FAILURE);
    }
  
  if(!read_factors(factors, q, facs_file))
    {
      printf("Fatal error reading facs file. Exiting\n");
      exit(FAILURE);
    }
  
  outfile=fopen(argv[7],"wb");
  if(!outfile)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[5]);
      exit(1);
    }

  if(!(Fft_vec=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*N)))
    {
      printf("Fatal error allocating memory for Fft_vec. Exiting.\n");
      exit(FAILURE);
    }

  if(!(ws=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*N/2)))
    {
      printf("Fatal error allocating memory for ws. Exiting.\n");
      exit(FAILURE);
    }

  initfft(N,ws);
  
  for(i=0;i<N;i++)
    mpfi_c_init(Fft_vec[i]);
  /*
  for(i=0;i<N;i++)
    mpfi_c_set_ui(Fft_vec[i],i,0);
  nifft(Fft_vec,N,ws);

  for(i=0;i<N;i++)
    mpfi_c_print(Fft_vec[i]);
  exit(0);
  */

  printf("Computing g_e(m).\n");
  g_e(t0,h,B,N,20,Fft_vec);
  for(i=0;i<N;i++) {printf("g_e(%d)=",i); mpfi_c_printn(Fft_vec[i],20);}
  /*
  printf("Computing g_e_twiddle(m).\n");
  g_e_twiddle(h,B,N,Fft_vec);
  for(i=0;i<N;i++) {printf("g_e_twiddle(%d)=",i); mpfi_c_printn(Fft_vec[i],20);}
  */
  printf("Computing G_e_twiddle(n).\n");
  G_e_twiddle(N,ws,Fft_vec);
  for(i=0;i<N;i++) {printf("G_e_twiddle(%d)=",i); mpfi_c_printn(Fft_vec[i],20);}
  exit(0);

  printf("Computing G_e(n).\n");
  G_e(A,h,t0,N,Fft_vec);
  for(i=0;i<N;i++) {printf("G_e(%d)=",i); mpfi_c_printn(Fft_vec[i],10);}
  

  exps=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*q);
  if(!exps)
    {
      printf("Fatal error allocating memory for exps. Exiting.\n");
      exit(FAILURE);
    }
  for(i=1;i<q;i++)
    mpfi_c_init(exps[i]);
  calc_exps(q,exps);
  if(!(chis=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*q)))
    {
      printf("Failed to allocate memory for chis.Exiting\n");
      exit(FAILURE);
    }

  for(i=0;i<q;i++)
    mpfi_c_init(chis[i]);

  mpfi_c_init(omega);
  mpfi_init(sqrt_q);
  mpfi_set_ui(sqrt_q,q);
  mpfi_sqrt(sqrt_q,sqrt_q);
  for(i=0;i<factors[q].phi;i++)
    if(primitive(q,i,factors))      
      {
        make_chis(q,i,factors,chis);
        calc_omega(omega,q,chis,exps,sqrt_q);
      }

  fclose(outfile);
}



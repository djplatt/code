#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "../includes/mpfi.h"
#include "../includes/mpfi_io.h"
#include "../includes/mpfi_c.h"
//#include "./f_defs.h"
#define one_over_A ((double) 5.0/64.0)
#define LN_SMALL ((double) -230.0) // exp(LN_SMALL) is tiny but a normaised float 
//#define TAYLOR_TERMS (30)  // how many terms to use in Taylor Series approx
//#define MIN_M_FOR_FFT (50) // if no. f_hat terms required less than this, use simple sum
#define TARGET_ACC (100) // aiming for F_hat_err < exp(-TARGET_ACC) used 40 for main GRH run

//#include "./qt.h"

//#define EVEN
#ifdef EVEN
#define E_FAC ((double) 4.0)
#define F1 ((double) 0.5)
#define F2 ((double) 0.25) 
#else
#define E_FAC ((double) 7.03)
#define F1 ((double) 1.5)
#define F2 ((double) 0.75) 
#endif

#define my_fwrite(a,b,c,d) {if(fwrite(a,b,c,d)!=c) { printf("Error writing data to file. Exiting.\n");exit(0);}}
#define fatal_error(str) {printf("%s. Exiting.\n",str);exit(0);}

#define MAX_Q (1500000) // actually 750000 for odd q
#define MAX_CONV (1<<21)//23) // smallest pwr of 2 >=phi(Q)-1
// this is sufficient because for prime q we do two half length Bluesteins
#define MAX_FACS (8) // 2*3*5*7*11*13*17>MAX_Q
#define MAX_DIMS (MAX_FACS+1) // plus 1 for 2^n trick

//#define MAX_SIMPLE_DFT (10000) // use simple O(n^2) algorithm for n<= this
// defined in mpfi_fft.h

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


#include "../includes/mpfi_fft.h"

/*
#define MAX_Q (2*3*5*7*11*13*19)
#define MAX_FACS (7)

#define true (1==1)
#define false (1==0)

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
/*
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
*/

inline void X_x(mpfi_ptr res, mpfi_ptr x, mpfi_ptr delta, mpfi_ptr logdelta, mpfi_ptr logq)
{
  mpfi_sqr(res,x);
  mpfi_add(res,res,logdelta);
  mpfi_add(res,res,mpfi_ln_pi);
  mpfi_sub(res,res,logq);
  mpfi_sub(res,res,delta);
  mpfi_exp(res,res);
  //return(exp(x*2+logdelta+d_log_pi-logq-delta));
}

inline void log_X_x(mpfi_ptr res, mpfi_ptr x, mpfi_ptr delta, mpfi_ptr logdelta, mpfi_ptr logq)
{
  mpfi_sqr(res,x);
  mpfi_add(res,res,logdelta);
  mpfi_add(res,res,mpfi_ln_pi);
  mpfi_sub(res,res,logq);
  mpfi_sub(res,res,delta);
  //return(x*2+logdelta+d_log_pi-logq-delta);
}

/*
#ifdef EVEN
void E (mpfi_ptr res, mpfi_ptr t, mpfi_ptr eta_pi_by_4, mpfi_ptr logq)
{
  int_complex lng=lngamma(int_complex(int_double(0.25),t/2));
  int_double res=exp(lng.real+eta_pi_by_4*t+logq*0.3125-d_log_pi*0.5625+d_log_zeta_9_by_8+log(t*t+9.0/4.0)*0.15625);
  //print_int_double_str("E_o called with t=",t);
  //print_int_double_str("Returning",res);
  return(res);
}

int_double beta(const int_double &t)
{
  int_double res=d_pi/4-atan2(d_one,abs(t)*2)*0.5-d_one/d_pi/d_pi/abs(t*2-1.0/4.0)*4;
  return(res);
}

int_double F_twiddle_err(const unsigned int m, const int_double &eta_pi_by_4, const int_double &logq)
{
  int_double t1=int_double(m)/A+B;
  int_double t2=int_double(m)/A-B;
  int_double beta1=beta(t1);
  assert(beta1.left>-eta_pi_by_4.right);
  int_double beta2=beta(t2);
  assert(beta2.left>-eta_pi_by_4.right);
  int_double res=E(t1,eta_pi_by_4,logq)/(d_one-exp(-(beta1-eta_pi_by_4)*B));
  res+=E(t2,eta_pi_by_4,logq)/(d_one-exp(-(beta2+eta_pi_by_4)*B));
  res.left=res.right;
  return(res);
}

int_complex F_hat_twiddle_err( const int_double &delta, const int_double &logdelta, const int_double &logq)
{
  int_double X_Pi_A=X_x(pi_A,delta,logdelta,logq);
  assert(X_Pi_A.left>1.0);
  int_double log_res=log(int_double(8))+two_pi_A-X_Pi_A+log(int_double(0.5)/X_Pi_A+1.0)-logq*0.25-logdelta*0.5-log(exp_minus_pi_A);
  if(log_res.left<-LN_SMALL)
    return(c_small);
  log_res=exp(log_res);
  log_res.left=log_res.right;
  return(int_complex(log_res,log_res));
}

#else

int_double E (const int_double &t, const int_double &eta_pi_by_4, const int_double &logq)
{
  int_complex lng=lngamma(int_complex(int_double(0.75),t/2));
  int_double res=exp(lng.real+eta_pi_by_4*t+logq*0.3125-d_log_pi*1.0625+d_log_zeta_9_by_8+log(t*t+9.0/4.0)*0.15625);
  //print_int_double_str("E_o called with t=",t);
  //print_int_double_str("Returning",res);
  return(res);
}

int_double beta(const int_double &t)
{
  int_double res=d_pi/4-atan2(d_one,abs(t)*2)*1.5-d_one/d_pi/d_pi/abs(t*2-9.0/4.0)*4;
  return(res);
}



int_double F_twiddle_err(const unsigned int m, const int_double &eta_pi_by_4, const int_double &logq)
{
  int_double t1=int_double(m)/A+B;
  int_double t2=int_double(m)/A-B;
  int_double beta1=beta(t1);
  assert(beta1.left>-eta_pi_by_4.right);
  int_double beta2=beta(t2);
  assert(beta2.left>-eta_pi_by_4.right);
  int_double res=E(t1,eta_pi_by_4,logq)/(d_one-exp(-(beta1-eta_pi_by_4)*B));
  res+=E(t2,eta_pi_by_4,logq)/(d_one-exp(-(beta2+eta_pi_by_4)*B));
  res.left=res.right;
  return(res);
}




int_complex F_hat_twiddle_err( const int_double &delta, const int_double &logdelta, const int_double &logq)
{
  int_double X_Pi_A=X_x(pi_A,delta,logdelta,logq);
  assert(X_Pi_A.left>1.0);
  int_double log_res=log(int_double(8))+two_pi_A*3.0-X_Pi_A+log(int_double(0.5)/X_Pi_A+1.0)*1.5-logq*0.75-logdelta*0.5-log(exp_minus_pi_A);
  if(log_res.right>-LN_SMALL)
    return(c_small);
  log_res=exp(log_res);
  log_res.left=log_res.right;
  return(int_complex(log_res,log_res));
}
#endif
*/

inline long unsigned int phi(long unsigned int p, long unsigned int p_n)
{
  return(p_n-p_n/p);
}

long unsigned int pow_mod(long unsigned int a, long unsigned int b, long unsigned int m)
{
  long unsigned int a_pw=a,pw=b,res=1;
  while(true)
    {
      //printf("pw=%ld a_pw=%ld res=%ld\n",pw,a_pw,res);
      if(pw&1)
	res=(res*a_pw)%m;
      pw>>=1;
      if(pw==0)
	return(res);
      a_pw=(a_pw*a_pw)%m;
    }
}
      
unsigned long int pr(unsigned long int i, factor *factors)
{
  long unsigned int phi=factors[i].phi,p,good,j;
  for(p=2;p<i;p++)
    {
      if(gcd(p,i)!=1)
	continue;
      good=true;
      for(j=0;j<factors[phi].num_facs;j++)
	{
	  if(pow_mod(p,phi/factors[phi].primes[j],i)!=1)
	    continue;
	  good=false;
	  break;
	}
      if(good)
	return(p);
    }
}

int make_factors(factor *factors, long unsigned int q_end)
{
  printf("Making factor database upto %lu.\n",q_end);
  long unsigned int *primes,max_p=floor(sqrt(q_end)),i,j,f,n;
  primes=(long unsigned int *) malloc(sizeof(long unsigned int)*(q_end+1));
  for(i=2;i<=q_end;i++)
    primes[i]=0;
  for(i=4;i<=q_end;i+=2)
    primes[i]=2;
  for(i=3;i<=max_p;i+=2)
    if(primes[i]==0)
      for(j=i*i;j<=q_end;j+=i)
	if(primes[j]==0)
	  primes[j]=i;
  printf("Prime sieve completed.\n");
  // now each entry primes[i] is 0 if i prime, else = smallest prime factor
  factors[3].num_facs=1;
  factors[3].primes[0]=3;
  factors[3].facs[0]=3;
  factors[4].num_facs=1;
  factors[4].primes[0]=2;
  factors[4].facs[0]=4;

  for(f=5;f<=q_end;f++)
    if(primes[f]==0) // a prime
      {
	factors[f].num_facs=1;
	factors[f].primes[0]=f;
	factors[f].facs[0]=f;
      }
    else
      {
	factors[f].primes[0]=primes[f];
	n=f/primes[f];
	if(factors[n].primes[0]==primes[f])
	  {
	    factors[f].num_facs=factors[n].num_facs;
	    factors[f].facs[0]=primes[f]*factors[n].facs[0];
	    for(i=1;i<factors[n].num_facs;i++)
	      {
		factors[f].primes[i]=factors[n].primes[i];
		factors[f].facs[i]=factors[n].facs[i];
	      }
	  }
	else
	  {
	    factors[f].num_facs=factors[n].num_facs+1;
	    factors[f].facs[0]=primes[f];
	    factors[f].primes[0]=primes[f];
	    for(i=1;i<factors[f].num_facs;i++)
	      {
		factors[f].primes[i]=factors[n].primes[i-1];
		factors[f].facs[i]=factors[n].facs[i-1];
	      }
	  }
      }

  free(primes);
  printf("Factors computed.\n");
  // now calculate phi(f)
  for(i=3;i<=q_end;i++)
    {
      factors[i].phi=1;
      for(j=0;j<factors[i].num_facs;j++)
	factors[i].phi*=phi(factors[i].primes[j],factors[i].facs[j]);
    }
  printf("phi computed.\n");

  //now do the prim roots
  factors[3].pr=2;
  factors[4].pr=3;
  //long unsigned int max_pr=3;
  for(i=5;i<=q_end;i++)
    {
      if(((factors[i].num_facs==1)&&(factors[i].primes[0]!=2))|| // p^m, p an odd prime
	 ((factors[i].num_facs==2)&&(factors[i].facs[0]==2)))    // 2p^m
	{
	  factors[i].pr=pr(i,factors);
	  /*
	  if(factors[i].pr>max_pr)
	    {
	      max_pr=factors[i].pr;
	      printf("New largest pr=%lu for modulus %ld\n",max_pr,i);
	    }
	  */
	}
      else
	factors[i].pr=0;
    }
  printf("pr's computed.\n");
  /*  
  for(int i=3;i<50;i++)
    {
      printf("%ld has %ld factors, phi(%ld)=%ld, pr(%ld)=%ld, factors are",i,factors[i].num_facs,i,factors[i].phi,i,factors[i].pr);
      for(int j=0;j<factors[i].num_facs;j++)
	printf(" %ld %ld",factors[i].primes[j],factors[i].facs[j]);
      printf("\n");
    }
  */
  return(true);
}


void print_usage(char *str)
{
  printf("Usage:- %s prec q H.\n",str);
  exit(1);
}

double QT;

unsigned long int calc_N(unsigned long int q)
{
  unsigned long int res=1;
  double target=(ceil(QT/(10.0*q))*10.0+30.0)/one_over_A*7.5;
  while(res<target) res<<=1;
  return(res);
}

// results in attenuation of about 1/100 at height QT/q+30
double calc_eta(unsigned int q)
{
  return(1-E_FAC/(QT/q+30.0)); // odd = 7.03
}

//#define TARGET_ACC (100) // aiming for F_hat_err < exp(-TARGET_ACC)

inline unsigned int calc_M(double x, double log2byq, double cosetaetc)
{
  double b=log2byq+x*F1; // odd =1.5
  //printf("x=%e\nlog(2)-log(q)*0.75=%e\ncos(eta)=%e\n",x,log2byq,cosetaetc);
  //printf("b=%e\n",b);
  //printf("lambda=%e\n",exp(2*x)*cosetaetc);
  if(b>TARGET_ACC)
    return(1);
  return((unsigned int) sqrt((TARGET_ACC-b)/exp(2*x)/cosetaetc)+1);
}

mpfi_t b,a,r,lambda,d_small;
mpfi_c_t inner_exp,c_small;
#ifndef EVEN
#define MAX_LOG_N (1000)
mpfi_t log_n[MAX_LOG_N+1];
#endif

inline void F_hat_err(mpfi_ptr res, mpfi_ptr x, mpfi_ptr pi_sindelta_by_q, mpfi_ptr mlog2byq, unsigned int m)
{
  mpfi_add(lambda,x,x);
  mpfi_exp(lambda,lambda);
  mpfi_mul(lambda,lambda,pi_sindelta_by_q);
  //printf("lambda=");mpfi_print(lambda);
  mpfi_mul_d(b,x,F1);
  mpfi_add(b,b,mlog2byq);
  //printf("F1*x+log2-F2logq=");mpfi_print(b);
  mpfi_mul_ui(a,lambda,(m+1)*(m+1));
  //printf("(m+1)^2*lambda=");mpfi_print(a);
  mpfi_sub(a,b,a);
  mpfi_set_ui(b,m+1);
  mpfi_log(b,b);
  mpfi_add(a,a,b);
  //mpfi_exp(b,a);printf("a=");mpfi_print(b);
  mpfi_mul_d(r,lambda,(-3.0-2.0*m));
  mpfi_exp(r,r);
  // this is the n that only appears in odd characters
#ifndef EVEN
  mpfi_set_ui(res,2+m);
  mpfi_div_ui(res,res,1+m);
  mpfi_mul(r,r,res);
#endif
  mpfi_sub_d(r,r,1.0);
  mpfi_neg(r,r);
  //printf("r=");mpfi_print(r);
  mpfi_log(r,r);
  mpfi_sub(res,a,r);
  /*
  if(mpfi_cmp_d(res,LN_SMALL)<0)
    mpfi_set(res,d_small);
  else
    {
  */
      mpfi_exp(res,res);
      mpfi_neg(a,res);
      mpfi_put(res,a);
      //}
}

void F_hat_term(mpfi_c_ptr res, mpfi_c_ptr outer_exp, mpfi_c_ptr exp_2_u_x, unsigned int n)
{
  //printf("F_hat_o_term called with\n");
  //printf("outer_exp=");mpfi_c_print(outer_exp);
  //printf("exp_2_u_x=");mpfi_c_print(exp_2_u_x);
  //exit(0);
  mpfi_c_mul_ui(inner_exp,exp_2_u_x,n);
  mpfi_c_mul_ui(inner_exp,inner_exp,n);
  mpfi_c_sub(inner_exp,outer_exp,inner_exp);
#ifndef EVEN
  if(n>MAX_LOG_N)
    {
      mpfi_set_ui(log_n[0],n);
      mpfi_log(log_n[0],log_n[0]);
      mpfi_add(inner_exp->re,inner_exp->re,log_n[0]);
    }
  else
    mpfi_add(inner_exp->re,inner_exp->re,log_n[n]);
#endif
  if(mpfi_cmp_d(inner_exp->re,LN_SMALL)<0)
    mpfi_c_set(res,c_small);
  else
    mpfi_c_exp(res,inner_exp);
  //printf("F_hat_term Returning- ");mpfi_c_print(res);
}

inline unsigned int min(unsigned int x, unsigned int y)
{
  if(x<y) return(x);
  return(y);
}

void fill_q_ns (long unsigned int *q_n,// unsigned int *q_neg_n,
		long unsigned int pr, long unsigned int phi_q,
		long unsigned int q)
{
  long unsigned int i,pr_n;
  q_n[1]=0;
  q_n[pr]=1;
  pr_n=pr;
  for(i=2;i<phi_q;i++)
    {
      pr_n=(pr_n*pr)%q;
      q_n[pr_n]=i;
    }
}


void fill_q_ns_2s(long unsigned int *q_n, long unsigned int q)
{
  long unsigned int i,pr;
  pr=1;
  for(i=0;i<(q>>2);i++)
    {
      q_n[pr]=i;
      q_n[q-pr]=i+(q>>2);
      pr=(pr*5)%q;
    }
}

// on exit, offset[i] is where the i'th fraction a/q goes in the fft vector
void make_offsets(long unsigned int q, long unsigned int *q_n, long unsigned int *offset_vec, factor *factors)
{
  unsigned long int i,j,fac,fac1,offset,offset1,ptr;
  unsigned long int *a_n=(unsigned long int *)malloc(sizeof(long unsigned int)*factors[q].phi);
  if(!a_n)
    fatal_error("Failed to allocate memory for a_n.");

  for(i=1,j=0;i<q;i++)
    if(co_prime(i,q))
      a_n[j++]=i;
  
  ptr=0;   // ptr into offset vec
  fac=factors[q].facs[0];
  for(i=0;i<factors[q].phi;i++)
    {
      offset=q_n[a_n[i]%fac];
      offset1=fac;
      j=1;
      while(j<factors[q].num_facs)
	{
	  fac1=factors[q].facs[j];
	  offset*=factors[fac1].phi;
	  offset+=q_n[a_n[i]%fac1+offset1];
	  offset1+=fac1;
	  j++;
	}
      offset_vec[ptr++]=offset;
    }
  free(a_n);
}

long unsigned int make_lookup(long unsigned int *lookup, long unsigned int q, 
			      factor *factors, long unsigned int *dims,
			      long unsigned int *lookup1)
{
  long unsigned int phi=factors[q].phi,no_dims,fac,i,offset;
  int power_2;
  long unsigned int pr=factors[q].pr;
  if(pr==0) // no primitive root, bugger
    {
      no_dims=factors[q].num_facs;
      fac=factors[q].facs[0];        // the first p^n
      power_2=(factors[fac].pr==0);  // no conductor => p=2, n>=3
      if(power_2)
	{
	  no_dims++;
	  fill_q_ns_2s(lookup1,fac);    // use the {-1,1}X{5} trick
	  for(i=1;i<factors[q].num_facs;i++) // move all the dimensions up one
	    dims[i+1]=factors[factors[q].facs[i]].phi;
	  dims[1]=factors[q].facs[0]/4;
	  dims[0]=2;                         // slip in a two
	}
      else
	{
	  fill_q_ns(lookup1,factors[fac].pr,factors[fac].phi,fac); // use the generator
	  for(i=0;i<factors[q].num_facs;i++)
	    dims[i]=factors[factors[q].facs[i]].phi;
	}
      offset=fac;
      for(i=1;i<factors[q].num_facs;i++)  // do the rest on the factors, all will have generators
	{
	  fac=factors[q].facs[i];      
	  pr=factors[fac].pr;
	  fill_q_ns(lookup1+offset,pr,factors[fac].phi,fac);  // use the generator
	  offset+=fac;
	}
      make_offsets(q,lookup1,lookup,factors);    // reverse q_n so we know how to populate fft vector
      return(no_dims);
    }
  else
    {
      dims[0]=phi;
      fill_q_ns(lookup1,pr,phi,q);
      /*
      for(i=0;i<phi;i++)
	{
	  printf("lookup1[%lu]=%lu\n",i,lookup1[i]);
	  lookup[lookup1[i]]=i;
	}
      */
      //printf("calling make_offsets.\n");
      make_offsets(q,lookup1,lookup,factors);
      //printf("returned from make_offsets.\n");
      /*
  lookup[0]=0;
  long unsigned int j;
  for(i=1,j=0;i<q;i++)
    if(co_prime(i,q))
      lookup1[i]=j++;
  for(j=1,i=pr;j<phi;j++,i=(i*pr)%q)
    {
      //printf("j=%lu i=%lu lookup1[i]=%lu\n",j,i,lookup1[i]);
      lookup[j]=lookup1[i];
    }
  */
      return(1); // one dimensional DFT
    }
}

int primp(unsigned long int *coords, factor *factors)
{
  unsigned long int j;
  if((factors[0].primes[0]==2)&&(factors[0].facs[0]>4)) // 2^alpha with alpha>=3
    {
      /*
      if(coords[0]==0) // the first length 2 dim
	return(false); // is irrelevant
      */
      if((coords[1]&1)==0) // the second length 2^(alpha-2) dim
	return(false);
      for(j=2;j<factors[0].num_facs+1;j++) // the remaining dims
	if((coords[j]%factors[0].primes[j-1])==0)
	  return(false);
      return(true);
    }
  else
    {
      for(j=0;j<factors[0].num_facs;j++)
	if(coords[j]%factors[0].primes[j]==0)
	  return(false);
      return(true);
    }
}

int negp(unsigned long int *coords, factor *factors, unsigned long int n_dims)
{
  
  unsigned long int n,sum=0;
  for(n=0;n<n_dims;n++)
    sum+=coords[n];

  if(n_dims!=factors[0].num_facs) // 2^alpha with alpha>=3
    return((sum&1)==0);
  else
    return((sum&1)==1);
}
 

void make_err(mpfi_c_ptr res, mpfi_ptr err)
{
  mpfi_neg(res->re,err);
  mpfi_put(res->re,err);
  mpfi_set(res->im,res->re);
}

int main(int argc, char **argv)
{
  long int q1,prec;
  long unsigned int q,M,n,i,j,k,N;
  double eta,log2byq,cosetaetc,B;
  mpfi_t x,delta,logdelta,pi_sindelta_by_q,err,log2,logq;
  mpfi_t two_pi_by_B,eta_pi_by_4,mlog2byq,logpi,pibyq,pi;
  mpfi_c_t *res,term,outer_exp,exp_2_u_x,u_x,cerr;
  factor *factors;

  printf("Command Line:-");
  for(q=0;q<argc;q++)
    printf(" %s",argv[q]);
  printf("\n");

  if(argc!=4)
    print_usage(argv[0]);

  prec=atoi(argv[1]);
  if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    print_usage(argv[0]);
  
  mpfi_c_setup(prec);

  q1=atol(argv[2]);
  if(q1<3)
    print_usage(argv[0]);
  //printf("q=%d\n",q1);
  q=(long unsigned int) q1;
  if(q>MAX_Q)
    fatal_error("q exceeded max q.");
  if((q&3)==2)
    fatal_error("there are no primitive characters mod that q.");
  if(q1<6) q1=6; else q1=q+1;
  factors=(factor *)malloc(sizeof(factor)*q1);
  if(!factors)
    fatal_error("Failed to allocate memory for factors.");
  if(!make_factors(factors,q1-1))
    fatal_error("Error building factors database.");
  unsigned long int phi_q=factors[q].phi;
  //printf("phi(%lu)=%lu.\n",q,phi_q);

  unsigned long int dims[MAX_DIMS],dim;
  unsigned long int coords[MAX_DIMS];

  unsigned long int *lookup,*lookup1;
  lookup=(unsigned long int *)malloc(sizeof(long unsigned int)*phi_q);
  if(!lookup)
    fatal_error("Failed to allocate memory for lookup.");
  lookup1=(unsigned long int *)malloc(sizeof(long unsigned int)*q);
  if(!lookup1)
    fatal_error("Failed to allocate memory for lookup1.");
  dim=make_lookup(lookup,q,factors,dims,lookup1);
  free(lookup1);
  //for(i=0;i<phi_q;i++) printf("lookup[%lu]=%lu\n",i,lookup[i]);

  //printf("Returned from lookup.\n");
  
  for(i=0;i<dim;i++)
    coords[i]=0;

  int *isneg,*isprim;
  isneg=(int *)malloc(sizeof(int)*phi_q);
  isprim=(int *)malloc(sizeof(int)*phi_q);
  if(!(isneg&&isprim))
    fatal_error("Failed to allocate memory for isneg and/or isprim.");
  for(i=0;i<phi_q;i++)
    {
      isprim[i]=primp(coords,&factors[q]);
      isneg[i]=negp(coords,&factors[q],dim);
      unsigned long int last_dim=dim-1;
      coords[last_dim]++;
      while(coords[last_dim]==dims[last_dim])
	{
	  coords[last_dim]=0;
	  if(last_dim>0)
	    {
	      last_dim--;
	      coords[last_dim]++;
	    }
	  else
	    break;
	}
    }
  //for(i=0;i<phi_q;i++) if(isprim[i]) {if(isneg[i]) printf("%lu is negative.\n",i); else printf("%lu is positive.\n",i);}
  mpfi_c_t *omegas=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*phi_q);
  if(!omegas)
    fatal_error("Failed to allocate memory for omegas.");
  for(i=0;i<phi_q;i++)
    mpfi_c_init(omegas[i]);
  //printf("omegas initialised.\n");
  simple_dft_init();
  //printf("dft initialised calling prep omegas with q=%lu.\n",q);
  prep_omegas_nd(lookup,q,omegas,dim,dims,phi_q);
  
  
  long unsigned int num_chars=0;
  for(i=0;i<phi_q;i++)
#ifdef EVEN
    if(isprim[i]&&(!isneg[i]))
#else
      if(isprim[i]&&isneg[i])
#endif
      {
	//printf("Calling finish omega on %lu\n",i);
	//printf("prepped omegas[%lu]=",i);mpfi_c_print(omegas[i]);
	finish_omega(omegas[i],isneg[i]);
	//if(isneg[i])
	//printf("neg ");
	//else
	//printf("pos ");
	//printf("finished omegas[%lu]=",i);mpfi_c_print(omegas[i]);
	num_chars++;
      }
  printf("There are %lu primitive ",num_chars);
#ifdef EVEN
  printf("even");
#else
  printf("odd");
#endif
  printf(" characters to consider.\n");

  QT=q*atof(argv[3]);
  eta=calc_eta(q);
  printf("eta set to %e\n",eta);

  N=calc_N(q);
  printf("N=%d\n",N);
  mpfi_c_t *Fft_vec=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*N*num_chars);
  if(!Fft_vec)
    fatal_error("Failed to allocate memory for Fft_vec.");
  for(i=0;i<N*num_chars;i++)
    mpfi_c_init(Fft_vec[i]);
  mpfi_c_t *ws=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*N/2);
  if(!ws)
    fatal_error("Failed to allocate memory for ws.");
  B=N*one_over_A;
  mpfi_c_init(cerr);
  mpfi_init(a);
  mpfi_init(r);
  mpfi_init(lambda);
  mpfi_init(b);
  mpfi_init(x);
  mpfi_init(two_pi_by_B);
  mpfi_init(delta);
  mpfi_init(pi);
  mpfi_init(eta_pi_by_4);
  mpfi_const_pi(pi);
  mpfi_div_d(two_pi_by_B,pi,B/2.0);
  //write_mpfi(outfile,two_pi_by_B);
  mpfi_init(pibyq);
  mpfi_div_ui(pibyq,pi,q);
  mpfi_init(logpi);
  mpfi_log(logpi,pi);
  mpfi_mul_d(delta,pi,(1.0-fabs(eta))/2.0);
  //write_mpfi(outfile,delta);
  mpfi_init(logdelta);
  mpfi_log(logdelta,delta);
  //write_mpfi(outfile,logdelta);
  mpfi_mul_d(eta_pi_by_4,mpfi_pi_by_4,eta);
  //write_mpfi(outfile,eta_pi_by_4);
  mpfi_init(pi_sindelta_by_q);
  mpfi_sin(pi_sindelta_by_q,delta);
  mpfi_mul(pi_sindelta_by_q,pi_sindelta_by_q,pi);
  mpfi_div_ui(pi_sindelta_by_q,pi_sindelta_by_q,q);
  //write_mpfi(outfile,pi_sindelta_by_q);
  //printf("Pi*sin(delta)/q=");mpfi_print(pi_sindelta_by_q);
  mpfi_init(err);
  mpfi_c_init(term);
  mpfi_init(d_small);
  mpfi_set_d(d_small,LN_SMALL);
  mpfi_exp(d_small,d_small);
  mpfi_neg(err,d_small);
  mpfi_put(d_small,err);
  mpfi_c_init(c_small);
  mpfi_set(c_small->re,d_small);
  mpfi_set(c_small->im,d_small);
  res=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*phi_q);
  if(!res)
    fatal_error("Failed to allocate memory for res.");
  for(i=0;i<phi_q;i++)
    mpfi_c_init(res[i]);
  mpfi_c_init(outer_exp);
  mpfi_c_init(inner_exp);
  mpfi_c_init(u_x);
  mpfi_c_init(exp_2_u_x);
#ifndef EVEN
  mpfi_init(log_n[0]);
  mpfi_init(log_n[1]);
  mpfi_set_ui(log_n[1],0);
  for(n=1;n<=MAX_LOG_N;n++)
    {
      mpfi_init(log_n[n]);
      mpfi_set_ui(log_n[n],n);
      mpfi_log(log_n[n],log_n[n]);
    }
#endif
  mpfi_init(log2);
  mpfi_init(logq);
  mpfi_set_d(log2,2.0);
  mpfi_log(log2,log2);
  //write_mpfi(outfile,log2);
  mpfi_set_ui(logq,q);
  mpfi_log(logq,logq);
  //write_mpfi(outfile,logq);
  mpfi_mul_d(logq,logq,F2); // odd =0.75
  mpfi_init(mlog2byq);
  mpfi_sub(mlog2byq,log2,logq);
  cosetaetc=cos(M_PI*eta/2.0)*M_PI/q;
  log2byq=log((double) 2.0)-log((double) q)*F2; // odd =0.75
  mpfi_set(u_x->im,eta_pi_by_4);
  printf("Max M will be %lu\n",calc_M(0,log2byq,cosetaetc));
  for(n=0;n<=N/2;n++) // result of final DFT is real so use conjugacy
    {
      mpfi_mul_ui(x,two_pi_by_B,n);
      mpfi_set(u_x->re,x);
      mpfi_c_mul_d(exp_2_u_x,u_x,2.0);
      mpfi_c_exp(exp_2_u_x,exp_2_u_x);
      mpfi_c_mul_i(exp_2_u_x,exp_2_u_x,pibyq);
      //printf("exp(2U(x))*pi/q=");mpfi_c_print(exp_2_u_x);
      mpfi_c_mul_d(outer_exp,u_x,F1); // odd =1.5
      mpfi_add(outer_exp->re,outer_exp->re,log2);
      mpfi_sub(outer_exp->re,outer_exp->re,logq);
      M=calc_M(2.0*n*M_PI/B,log2byq,cosetaetc);
      //printf("M=%d\n",M);
      //my_fwrite(&M,sizeof(unsigned int),1,outfile);
      F_hat_err(err,x,pi_sindelta_by_q,mlog2byq,M);
      make_err(cerr,err);
      //write_mpfi(outfile,err);
      //printf("cerr=");mpfi_c_print(cerr);

      for(i=1,j=0;i<min(q,M+1);i++) // j counts characters
	if(co_prime(i,q))
	  {
	    mpfi_c_zero(res[lookup[j]]); // lookup[j] tells us where in dft to put it
	    for(k=i;k<=M;k+=q)
	      {
		F_hat_term(term,outer_exp,exp_2_u_x,k);
		mpfi_c_inc(res[lookup[j]],term);
	      }
	    //printf("res[%lu]=",lookup[j]);mpfi_c_print(res[lookup[j]]);
	    j++;
	    //write_mpfi_c(outfile,res);
	  }
      if(M==1)
	for(j=1;j<phi_q;j++)
	  mpfi_c_set(res[j],res[0]);
      else
	{
	  for(;i<q;i++)
	    if(co_prime(i,q))
	      mpfi_c_zero(res[lookup[j++]]);
	  ndft(res,phi_q,dim,dims);
	  //printf("ndft done.\n");
	}
      /*
      for(i=0;i<phi_q;i++)
	{
	  printf("res[%lu]=",i);
	  mpfi_c_print(res[i]);
	}
      */
      for(i=0,j=0;i<phi_q;i++)
#ifdef EVEN
	if(isprim[i]&&(!isneg[i]))
#else
	if(isprim[i]&&isneg[i])
#endif
	  {
	    /*
	    if(n==0)
	      {
		mpfi_c_init(omegas[i]);
		mpfi_c_conj(omegas[i],res[i]);
		mpfi_sqr(Fft_vec[j*N]->re,res[i]->re);
		mpfi_sqr(Fft_vec[j*N]->im,res[i]->im);
		mpfi_add(Fft_vec[j*N]->re,Fft_vec[j*N]->re,Fft_vec[j*N]->im);
		mpfi_sqrt(Fft_vec[j*N]->im,Fft_vec[j*N]->re);
		mpfi_div(omegas[i]->re,omegas[i]->re,Fft_vec[j*N]->im);
		mpfi_div(omegas[i]->im,omegas[i]->im,Fft_vec[j*N]->im);
	      }
	    */
	    mpfi_c_mul(Fft_vec[j*N+n],res[i],omegas[i]); // put this somewhere
	    mpfi_c_add(Fft_vec[j*N+n],Fft_vec[j*N+n],cerr); // add the error truncating sum at M
	    //printf("Fft[%lu][%lu]=",i,n);mpfi_c_print(Fft_vec[j*N+n]);
	    j++;
	  }
		  
    }
  // now we have num_chars length N FFT's to do
  // Fft_vec contains f_hat(n)
  // add an error to approximate f_hat->f_hat_twiddle
  // do the FFT
  // add an error to approximate f_twiddle->f
  printf("Going to perform FFT from f_hat_twiddle to f_twiddle.\n");
  initfft(N,ws);
  for(i=0;i<num_chars;i++)
    {
      for(n=1;n<N/2;n++)
	{
	  mpfi_set(Fft_vec[i*N+N-n]->re,Fft_vec[i*N+n]->re);
	  mpfi_neg(Fft_vec[i*N+N-n]->im,Fft_vec[i*N+n]->im);
	}
      // add f_hat_twiddle_error
      fft(Fft_vec+i*N,N,ws);
      if(!mpfi_has_zero(Fft_vec[N*i]->im))
	{
	  printf("Character %lu wasn't real.\n",i);
	  mpfi_c_print(Fft_vec[N*i]);
	}
      // add f_twiddle_error
    }
      
}


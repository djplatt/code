#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../includes/mpfi_c.h"
#include "../includes/mpfi_fft.h"
#include "../includes/fft_defs.h"
#include "../low_q/f_defs.h"

#define TRUE (0==0)
#define FALSE (1==0)
#define FAILURE (1)
#define debug printf("Reached line number %d\n",__LINE__)
#define bool int

unsigned int calc_N(unsigned int q)
{
  unsigned int res=1;
  double target=(ceil(QT/(10.0*q))*10.0+30.0)/one_over_A*5.0; // odd was *2.5
  while(res<target) res<<=1;
  return(res);
}

// results in attenuation of about 1/100 at height QT/q+30
double calc_eta(unsigned int q)
{
  return(1-4.0/(QT/q+30.0)); // odd = 7.03
}


// on entry s_vec[i] contains log((i+1)/sqrt(q))
// on exit                    "" - u_m
int calc_buck(int fft_len, int *buck, double B, mpfi_t *s_vec, int M, mpfi_ptr sqrt_q)
{
  int i,offset;
  double x;
  mpfi_t tmp;
  mpfi_init(tmp);
  //printf("s_vec[0]=");mpfi_printn(s_vec[0],10);
  mpfi_init(s_vec[0]);
  mpfi_set_ui(s_vec[0],1);
  mpfi_div(s_vec[0],s_vec[0],sqrt_q);
  mpfi_log(s_vec[0],s_vec[0]);
  x=mpfi_get_d(s_vec[0]);
  offset=mpfi_get_d(s_vec[0])*B-0.5;
  offset=-offset;
  buck[0]=0;

  for(i=0;i<M;i++)
    {
      mpfi_init(s_vec[i]);
      mpfi_set_ui(s_vec[i],i+1);
      mpfi_div(s_vec[i],s_vec[i],sqrt_q);
      mpfi_log(s_vec[i],s_vec[i]);
      x=mpfi_get_d(s_vec[i]);
      if(x>=0.0)
	buck[i]=x*B-0.5+offset;
      else
	buck[i]=x*B+0.5+offset;
      mpfi_set_ui(tmp,buck[i]);
      mpfi_div_d(tmp,tmp,B);
      mpfi_sub(s_vec[i],s_vec[i],tmp);
    }
  /*
  for(i=0;i<30;i++)
    printf("buck[%d]=%d\n",i,buck[i]);
  */
  assert(buck[M-1]<fft_len/2);
  for(i=1;i<M;i++)
    buck[i]=fft_len-buck[i];
  /*
  for(i=0;i<30;i++)
    printf("buck[%d]=%d\n",i,buck[i]);
  exit(0);
  */
  mpfi_clear(tmp);
  return(offset);
}


void print_usage()
{
  printf("Usage:- f_hat_even prec q this_file num_files facs_file outfile.\n");
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
mpfi_t theta,theta1;
mpfi_c_t chi,chi_neg;

void make_chis(unsigned int q, unsigned int index, factor * factors, mpfi_c_t *chis)
{
  unsigned int fac_ptr,fac,phi,pr,prn,i,pw;

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
}

inline unsigned int calc_M(double x, double log2byq, double cosetaetc)
{
  double b=log2byq+x*1.5;
  //printf("x=%e\nlog(2)-log(q)*0.75=%e\ncos(eta)=%e\n",x,log2byq,cosetaetc);
  //printf("b=%e\n",b);
  //printf("lambda=%e\n",exp(2*x)*cosetaetc);
  if(b>TARGET_ACC)
    return(1);
  return((unsigned int) sqrt((TARGET_ACC-b)/exp(2*x)/cosetaetc)+1);
}

mpfi_t **polyf;

void calc_polyf(unsigned int K)
{
  mpfi_t mtmp,fac;
  mpfi_init(mtmp);
  mpfi_init(fac);
  unsigned int k,coeff,k1_ptr=0,k2_ptr=2;
  if(!(polyf=(mpfi_t **) malloc(sizeof(mpfi_t *)*(K-1))))
    {
      printf("Failed to allocate memory for polyf. Exiting.\n");
      exit(1);
    }
  for(k=0;k<K-1;k++)
    {
      if(!(polyf[k]=(mpfi_t *) malloc(sizeof(mpfi_t)*(k+2))))
	{
	  printf("Failed to allocate memory for polyf. Exiting.\n");
	  exit(1);
	}
      for(coeff=0;coeff<k+2;coeff++)
	mpfi_init(polyf[k][coeff]);
    }

  mpfi_set_d(polyf[0][0],0.5);
  mpfi_set_d(polyf[0][1],2.0); // (2f+1/2)
  mpfi_set_ui(fac,1);
  for(k=1;k<K-1;k++)
    {
      mpfi_mul_d(polyf[k][0],polyf[k-1][0],0.5);
      for(coeff=1;coeff<=k;coeff++)
	{
	  mpfi_mul_d(polyf[k][coeff],polyf[k-1][coeff],coeff*2.0+0.5);
	  mpfi_mul_ui(mtmp,polyf[k-1][coeff-1],2);
	  mpfi_add(polyf[k][coeff],polyf[k][coeff],mtmp);
	}
      mpfi_mul_ui(polyf[k][k+1],polyf[k-1][k],2);
      mpfi_mul_ui(fac,fac,k+1);
      for(coeff=0;coeff<=k+1;coeff++)
	mpfi_div(polyf[k][coeff],polyf[k][coeff],fac);
    }
  mpfi_clear(fac);
  mpfi_clear(mtmp);
}

inline unsigned int min(unsigned int x, unsigned int y)
{
  if(x<y) return(x);
  return(y);
}

void test_chis(unsigned int q, factor *factors, mpfi_c_t *chis)
{
  unsigned int p1,p2;
  mpfi_c_t chi_prod;
  mpfi_c_init(chi_prod);
  for(p1=1;p1<q;p1++)
    if(co_primes[p1])
      for(p2=p1;p2<q;p2++)
	if(co_primes[p2])
	  {
	    mpfi_c_mul(chi_prod,chis[p1],chis[p2]);
	    mpfi_c_sub(chi_prod,chi_prod,chis[(p1*p2)%q]);
	    if(!mpfi_c_contains_zero(chi_prod))
	      {printf("bad chis - chi(%d).chi(%d)-chi(%d)=",p1,p2,(p1*p2)%q);mpfi_c_print(chi_prod);}
	  }
  mpfi_c_clear(chi_prod);
}

mpfi_c_t G_even_tmp,mpfi_c_small;

void G_even(mpfi_c_ptr res, mpfi_c_ptr u_ieta)
{
  mpfi_c_mul_ui(res,u_ieta,2);
  mpfi_c_div_ui(G_even_tmp,u_ieta,2);
  mpfi_c_exp(res,res);
  mpfi_c_mul_i(res,res,mpfi_pi);
  mpfi_c_sub(res,G_even_tmp,res);
  mpfi_add(res->re,res->re,mpfi_log_2);
  /*
  if(mpfi_cmp_d(res->re,LN_SMALL)<0)
    {
      mpfi_c_set(res,mpfi_c_small);
      return;
    }
  */
  mpfi_c_exp(res,res);
  return;
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

#define MAX_LOG_N (1000)

mpfi_t log_n[MAX_LOG_N+1];

inline int neg_one(unsigned int q, mpfi_c_t *chis)
{
  int res;
  mpfi_c_add_ui(chis[q-1],chis[q-1],1);
  res=mpfi_c_contains_zero(chis[q-1]);
  mpfi_c_sub_ui(chis[q-1],chis[q-1],1);
  return(res);
}

mpfi_c_t calc_s_temp;

inline void calc_s(mpfi_c_ptr res, mpfi_c_ptr chi, mpfi_ptr s_term, unsigned int pow)
{
  unsigned int i;
  //printf("in calc_s with\nchi=");mpfi_c_print(chi);printf("s_term=");mpfi_print(s_term);printf("pow=%d\n",pow);
  mpfi_zero(calc_s_temp->im);
  mpfi_set(calc_s_temp->re,s_term);
  for(i=1;i<pow;i++)
    mpfi_mul(calc_s_temp->re,calc_s_temp->re,s_term);
  mpfi_c_mul(calc_s_temp,calc_s_temp,chi);   // chi(n)/sqrt(n)
  //printf("adding ");mpfi_c_print(calc_s_temp);
  mpfi_c_add(res,res,calc_s_temp);
}

mpfi_c_t calc_f_term,calc_f_n;

inline void calc_f(mpfi_c_ptr res, mpfi_c_ptr G, mpfi_c_ptr f, mpfi_t *poly, unsigned int k)
{
  unsigned int i;
  //printf("in calc_f with G=");mpfi_c_print(G);
  //printf("with f=");mpfi_c_print(f);
  //printf("and with coefficients\n");
  //for(i=0;i<k+2;i++)
  //  {
  //    printf("%d ",i);mpfi_print(poly[i]);
  //  }
  mpfi_zero(res->im);
  mpfi_set(res->re,poly[0]);  // constant term
  mpfi_c_mul_i(calc_f_term,f,poly[1]);  // term in f
  mpfi_c_add(res,res,calc_f_term);
  mpfi_c_set(calc_f_n,f);
  for(i=2;i<k+2;i++)
    {
      mpfi_c_mul(calc_f_n,calc_f_n,f);
      mpfi_c_mul_i(calc_f_term,calc_f_n,poly[i]);
      mpfi_c_add(res,res,calc_f_term);
    }
  mpfi_c_mul(res,res,G);
  //printf("set res=");mpfi_c_print(res);
}


void FFT_f_hat_even(unsigned int q, double eta, unsigned int N,double B,unsigned int M,
		    factor *factors, mpfi_c_t *chis, mpfi_ptr sqrt_q, mpfi_c_t *exps, int n0, int n1)
{
  printf("In FFT_f_hat_even.\n");
  unsigned int n,i,j,k,index,*buck,prim1;
  int m,prim,real_p,offset,fft_len=(n1-n0)*2;
  double B_by_2pi=B/2.0/M_PI,d_sqrt_q=sqrt(q);
  mpfi_c_t *f_vec,*ws,*S_vec,*G_vec,*G0_vec,*res_vec,u_ieta,chi,omega;
  // ws contain omegas for FFT
  // S_vec contain S^k_m
  // G_vec contain G^k(u_m)/k!
  // G0_vec contains G(u_m)
  // res_vec contains sum over all k
  // f_vec contains -Pi*exp(4*n*Pi/B+Pi*eta*i/2)
  mpfi_t *s_term_vec,um,two_pi_by_B,*sqrt_n;
  // s_term_vec contains (log(n/sqrt(q))-u_m)
  mpfi_c_init(u_ieta);
  mpfi_c_init(chi);
  mpfi_init(um);
  mpfi_init(two_pi_by_B);
  mpfi_c_init(omega);
  mpfi_div_d(two_pi_by_B,mpfi_2_pi,B);
  printf("fft length %d.\n",fft_len);
  if(!(ws=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*fft_len/2)))
    {
      printf("Error allocating memory for ws. Exiting.\n");
      exit(1);
    }
  if(!(S_vec=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*fft_len)))
    {
      printf("Error allocating memory for S_vec. Exiting.\n");
      exit(1);
    }
  if(!(G_vec=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*fft_len)))
    {
      printf("Error allocating memory for G_vec. Exiting.\n");
      exit(1);
    }
  if(!(G0_vec=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*fft_len/2)))
    {
      printf("Error allocating memory for G0_vec. Exiting.\n");
      exit(1);
    }
  if(!(res_vec=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*fft_len/2)))
    {
      printf("Error allocating memory for res_vec. Exiting.\n");
      exit(1);
    }
  if(!(f_vec=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*fft_len/2)))
    {
      printf("Error allocating memory for f_vec. Exiting.\n");
      exit(1);
    }
  if(!(s_term_vec=(mpfi_t *) malloc(sizeof(mpfi_t)*M)))
    {
      printf("Error allocating memory for s_term_vec. Exiting.\n");
      exit(1);
    }
  if(!(buck=(unsigned int *) malloc(sizeof(mpfi_c_t)*M)))
    {
      printf("Error allocating memory for buck. Exiting.\n");
      exit(1);
    }
  if(!(sqrt_n=(mpfi_t *) malloc(sizeof(mpfi_t)*M)))
    {
      printf("Error allocating memory for sqrt_n. Exiting.\n");
      exit(1);
    }
  // sqrt[i]=sqrt(i+1)
  for(n=0;n<M;n++)
    {
      mpfi_init(sqrt_n[n]);
      mpfi_set_ui(sqrt_n[n],n+1);
      mpfi_sqrt(sqrt_n[n],sqrt_n[n]);
    }

  printf("Initialising FFT.\n");
  initfft(fft_len,ws);

  printf("Calculating buckets.\n");
  offset=calc_buck(fft_len,buck,B,s_term_vec,M,sqrt_q); // s_term_vec[i]=log((i+1)/sqrt(q))-u_m


  printf("Initialising vectors.\n");
  mpfi_mul_d(u_ieta->im,mpfi_pi_by_4,eta);
  for(i=0;i<fft_len/2;i++)
    {
      mpfi_mul_ui(u_ieta->re,two_pi_by_B,i);
      mpfi_c_init(f_vec[i]);
      mpfi_c_mul_ui(f_vec[i],u_ieta,2);
      mpfi_c_exp(f_vec[i],f_vec[i]);
      mpfi_c_mul_i(f_vec[i],f_vec[i],mpfi_pi);
      mpfi_c_neg(f_vec[i],f_vec[i]);
      mpfi_c_init(G0_vec[i]);
      G_even(G0_vec[i],u_ieta);
      mpfi_c_init(G_vec[i]);
      mpfi_c_init(S_vec[i]);
      mpfi_c_init(res_vec[i]);
    }
  for(;i<fft_len;i++)
    {
      mpfi_c_init(G_vec[i]);
      mpfi_c_init(S_vec[i]);
    }
  printf("Calculating polynomial in f.\n");
  calc_polyf(TAYLOR_TERMS); // calc G'(u) up to G^(TAYLOR_TERMS-1)(u)
  prim=-1;
  for(index=0;index<factors[q].phi;index++)
    {
      printf("Index = %d\n",index);
      if(!primitive(q,index,factors))
	continue;
      prim++;
      prim1=conj_j(prim,q,factors);
      if(prim1<prim)
	continue;
      real_p=(prim1==prim);
      printf("Making chis.\n");
      make_chis(q,index,factors,chis);
      printf("Testing chis.\n");
      test_chis(q,factors,chis); // warning O(q^2)
      // wasteful, should check negative or not before make_chis
      if(neg_one(q,chis))
	continue;
      calc_omega(omega,q,chis,exps,sqrt_q);
      for(i=0;i<offset;i++)
	mpfi_c_zero(G_vec[i]);
      for(i=0,j=offset;i<=fft_len/2;i++,j++)
	mpfi_c_set(G_vec[j],G0_vec[i]);
      for(;j<fft_len;j++)
	mpfi_c_zero(G_vec[j]);
      printf("Setting up S_vec.\n");
      for(i=0;i<fft_len;i++)
	mpfi_c_zero(S_vec[i]);
      for(n=1;n<=M;n++)
	if(co_primes[n%q])
	  {
	    mpfi_c_div_i(chi,chis[n%q],sqrt_n[n-1]);
	    mpfi_c_add(S_vec[buck[n-1]],S_vec[buck[n-1]],chi); // chi(n)/sqrt(n)(log(n/sqrt(q))-u_m)^0
	    //printf("S_vec[%d] set to",buck[n]);mpfi_c_print(S_vec[buck[n]]);
	  }
      printf("Doing G^(0) convolution.\n");
      convolve(G_vec,G_vec,S_vec,fft_len,ws); // G terms done
      printf("Done G^(0) convolution.\n");
      
      for(i=0;i<fft_len/2;i++)
	mpfi_c_set(res_vec[i],G_vec[i]);
      debug;

      for(k=0;k<TAYLOR_TERMS-1;k++)
	{
	  mpfi_c_set(G_vec[0],G0_vec[0]);
	  mpfi_c_zero(S_vec[0]);
	  for(i=1;i<fft_len/2;i++)
	    {
	      mpfi_c_zero(G_vec[i]);
	      mpfi_c_set(G_vec[fft_len-i],G0_vec[i]);
	      mpfi_c_zero(S_vec[i]);
	    }
	  //debug;
	  mpfi_c_zero(G_vec[fft_len/2]);
	  for(;i<fft_len;i++)
	    mpfi_c_zero(S_vec[i]);
	  
	  //debug;

	  /*
	  for(i=0;i<k+2;i++)
	    {printf("coeff[%d]=",i);mpfi_print(polyf[k][i]);}
	  */
	  calc_f(G_vec[0],G0_vec[0],f_vec[0],polyf[k],k);
	  //printf("G_vec[0]=");mpfi_c_print(G_vec[0]);
	  for(i=1;i<fft_len/2;i++)
	    {
	      calc_f(G_vec[fft_len-i],G0_vec[i],f_vec[i],polyf[k],k); // multiplies G(u) by P(f(u))
	      //printf("G_vec[%d]=",fft_len-i);mpfi_c_print(G_vec[fft_len-i]);
	    }
	  /*
	  for(i=fft_len/2;i<fft_len;i+=fft_len/32)
	    {
	      printf("G_vec[%d]=",i);
	      mpfi_c_print(G_vec[i]);
	    }
	  */
	  //debug;
	  for(n=1;n<=M;n++)
	    if(co_primes[n%q])
	      {
		mpfi_c_div_i(chi,chis[n%q],sqrt_n[n]);
		//debug;
		calc_s(S_vec[buck[n]],chi,s_term_vec[buck[n]],k+1);
		//debug;
	      }
	  printf("Doing G^(%d) convolution.\n",k+1);
	  convolve(G_vec,G_vec,S_vec,fft_len,ws);
	  printf("Done G^(%d) convolution.\n",k+1);
	  /*
	  for(n=1;n<fft_len/2;n<<=1)
	    {
	      printf("G_vec[%d]=",n);
	      mpfi_c_print(G_vec[n]);
	    }
	  */
	  for(n=0;n<fft_len/2;n++)
	    mpfi_c_add(res_vec[n],res_vec[n],G_vec[n]);
	  /*
	  for(n=0;n<fft_len/2;n+=fft_len/64)
	    {
	      printf("G^(%d) res[%d]=",k+1,n);mpfi_c_print(res_vec[n]);
	    }
	  */
	}
      for(n=0;n<fft_len/2;n++)
	mpfi_c_mul(res_vec[n],res_vec[n],omega);
      for(n=0;n<fft_len/2;n+=fft_len/64)
	{
	  printf("res[%d]=",n);
	  mpfi_c_print(res_vec[n]);
	}
      if(prim>5)
	exit(0);
    }
  mpfi_c_clear(u_ieta);
  mpfi_clear(um);
}

int main(int argc, char **argv)
{
  int q1,prec,this_file, num_files;
  unsigned int q,M,n,i,j,k,N,m;
  double eta,log2byq,cosetaetc,B;
  //mpfi_t x,delta,logdelta,pi_sindelta_by_q,err,log2,logq;
  //mpfi_t two_pi_by_B,eta_pi_by_4,mlog2byq,logpi,pibyq,pi;
  //mpfi_c_t res,term,outer_exp,exp_2_u_x,u_x;
  FILE *outfile,*facs_file;
  factor *factors;
  mpfi_c_t *chis,*exps;
  mpfi_t sqrt_q;


  if(argc!=7)
    print_usage();

  prec=atoi(argv[1]);
  if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    print_usage();
  printf("Running at %d bits of precision.\n",prec);
  //mpfr_set_default_prec(prec);
    
  mpfi_c_setup(prec);
  
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

  eta=calc_eta(q);
  printf("eta=%f\n",eta);
  N=calc_N(q);
  printf("N=%ld\n",N);
  B=N*one_over_A;
  printf("B=%f\n",B);
  this_file=atoi(argv[3]);
  num_files=atoi(argv[4]);
  printf("Processing file %d of %d\n",this_file,num_files);
  if((this_file<0)||(this_file>=num_files))
    print_usage();
  facs_file=fopen(argv[5],"r");
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
  
  outfile=fopen(argv[6],"wb");
  if(!outfile)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[5]);
      exit(1);
    }
  /*
  fwrite(&q,sizeof(unsigned int),1,outfile);
  fwrite(&eta,sizeof(double),1,outfile);
  fwrite(&n0,sizeof(unsigned int),1,outfile);
  fwrite(&n1,sizeof(unsigned int),1,outfile);
  fwrite(&N,sizeof(unsigned int),1,outfile);
  */
  if(!(chis=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*q)))
    {
      printf("Failed to allocate memory for chis.Exiting\n");
      exit(1);
    }
  for(i=0;i<q;i++)
    mpfi_c_init(chis[i]);
  if(!(exps=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*q)))
    {
      printf("Failed to allocate memory for exps.Exiting\n");
      exit(1);
    }
  for(i=1;i<q;i++)
    mpfi_c_init(exps[i]);
  calc_exps(q,exps);
/*
  for(i=1;i<q;i++)
      if(co_primes[i])
      {
	  printf("exps[%d]=",i);
	  mpfi_c_print(exps[i]);
      }
*/
  mpfi_init(theta);
  mpfi_init(theta1);
  mpfi_c_init(chi);
  mpfi_c_init(chi_neg);
  mpfi_c_init(mpfi_c_small);
  mpfi_c_init(G_even_tmp);
  mpfi_c_init(omega_term);
  mpfi_c_init(calc_s_temp);
  mpfi_c_init(calc_f_term);
  mpfi_c_init(calc_f_n);
  mpfi_set_d(mpfi_c_small->re,LN_SMALL);
  mpfi_exp(mpfi_c_small->re,mpfi_c_small->re);
  mpfi_neg(mpfi_c_small->im,mpfi_c_small->re);
  mpfi_put(mpfi_c_small->re,mpfi_c_small->im);
  mpfi_set(mpfi_c_small->im,mpfi_c_small->re);
  //printf("small=");mpfi_c_print(mpfi_c_small);
  mpfi_init(sqrt_q);
  mpfi_set_ui(sqrt_q,q);
  mpfi_sqrt(sqrt_q,sqrt_q);  

  mpfi_init(log_n[1]);
  mpfi_set_ui(log_n[1],0);
  for(n=2;n<=MAX_LOG_N;n++)
    {
      mpfi_init(log_n[n]);
      mpfi_set_ui(log_n[n],n);
      mpfi_div(log_n[n],log_n[n],sqrt_q);
      mpfi_log(log_n[n],log_n[n]);
    }

  M=calc_M(0.0,log((double) 2.0)-log((double) q)*0.25,cos(M_PI*eta/2.0)*M_PI/q);
  printf("We are using %d terms in f_hat sum.\n",M);
  FFT_f_hat_even(q,eta,N,B,M,factors,chis,sqrt_q,exps,N/2/num_files*this_file,N/2/num_files*(this_file+1));



  fclose(outfile);
  
}



#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../includes/mpfi_c.h"
#include "../includes/fft_defs.h"
#include "./f_defs.h"

#define TRUE (0==0)
#define FALSE (1==0)
#define FAILURE (1)
#define debug printf("Reached line number %d\n",__LINE__)
#define bool int

void print_usage()
{
  printf("Usage:- f_hat_even prec q eta n0 n1 N facs_file outfile.\n");
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

      //printf("theta=");mpfi_print(theta);
      if(pr==0)
	{
	  debug;
	  mpfi_div_ui(theta,mpfi_2_pi,phi/2);
	  mpfi_mul_ui(theta,theta,index);
	  mpfi_neg(theta,theta);
	  mpfi_print(theta);
	  prn=5;
	  for(pw=1;pw<phi/2;pw++)
	    {
	      mpfi_mul_ui(theta1,theta,pw);
	      mpfi_sin(chi->im,theta1);
	      mpfi_cos(chi->re,theta1);
	      for(i=prn;i<q;i+=fac)
		if(co_primes[i])
		  mpfi_c_mul(chis[i],chis[i],chi);
	      if(index&2)
		mpfi_c_neg(chi,chi);
	      for(i=fac-prn;i<q;i+=fac)
		if(co_primes[i])
		  mpfi_c_mul(chis[i],chis[i],chi);
	      prn=(prn*5)%fac;
	    }
	  if(index&2)
	    for(i=fac-1;i<q;i+=fac)
	      mpfi_c_neg(chis[i],chis[i]);
	}
      else
	{
	  mpfi_div_ui(theta,mpfi_2_pi,phi);
	  mpfi_mul_ui(theta,theta,index);
	  mpfi_neg(theta,theta);
	  prn=pr;
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
    }
}

inline unsigned int calc_M(double x, double log2byq, double cosetaetc)
{
  double b=log2byq+x*1.5;
  printf("x=%e\nlog(2)-log(q)*0.75=%e\ncos(eta)=%e\n",x,log2byq,cosetaetc);
  printf("b=%e\n",b);
  printf("lambda=%e\n",exp(2*x)*cosetaetc);
  if(b>TARGET_ACC)
    return(1);
  return((unsigned int) sqrt((TARGET_ACC-b)/exp(2*x)/cosetaetc)+1);
}
/*
mpfi_t b,a,r,lambda,d_small;
mpfi_c_t inner_exp,c_small;

inline void F_hat_err(mpfi_ptr res, mpfi_ptr x, mpfi_ptr pi_sindelta_by_q, mpfi_ptr mlog2byq, unsigned int m)
{
  mpfi_add(lambda,x,x);
  mpfi_exp(lambda,lambda);
  mpfi_mul(lambda,lambda,pi_sindelta_by_q);
  //printf("lambda=");mpfi_print(lambda);
  mpfi_mul_d(b,x,1.5);
  mpfi_add(b,b,mlog2byq);
  //printf("1.5x+log2-0.75logq=");mpfi_print(b);
  mpfi_mul_ui(a,lambda,(m+1)*(m+1));
  //printf("(m+1)^2*lambda=");mpfi_print(a);
  mpfi_sub(a,b,a);
  mpfi_set_ui(b,m+1);
  mpfi_log(b,b);
  mpfi_add(a,a,b);
  //mpfi_exp(b,a);printf("a=");mpfi_print(b);
  mpfi_mul_d(r,lambda,(-3.0-2.0*m));
  mpfi_exp(r,r);
  mpfi_mul_d(r,r,(2.0+m)/(1.0+m));
  mpfi_sub_d(r,r,1.0);
  mpfi_neg(r,r);
  //printf("r=");mpfi_print(r);
  mpfi_log(r,r);
  mpfi_sub(res,a,r);
  if(mpfi_cmp_d(res,LN_SMALL)<0)
    mpfi_set(res,d_small);
  else
    {
      mpfi_exp(res,res);
      mpfi_neg(a,res);
      mpfi_put(res,a);
    }
}
*/
mpfi_t **polyf;

void calc_polyf(unsigned int K)
{
  mpfi_t mtmp;
  mpfi_init(mtmp);
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
    }
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
  if(mpfi_cmp_d(res->re,LN_SMALL)<0)
    {
      mpfi_c_set(res,mpfi_c_small);
      return;
    }
  mpfi_c_exp(res,res);
  return;
}

unsigned int conj_j (unsigned int j, unsigned int num_chi, unsigned int q)
{
	if(q&7) // q not divisible by 8
	  return(num_chi-j-1);
	if(j<(num_chi>>1))
	  return((num_chi>>1)-j-1);
	return(num_chi+(num_chi>>1)-j-1);
}

bool primitive(unsigned int q, unsigned int index, factor *factors)
{
  int f;
  for(f=factors[q].num_facs-1;f>=0;f--)
    {
      if(((index%factors[factors[q].facs[f]].phi)%factors[q].primes[f])==0)
	return(FALSE);
      index/=factors[factors[q].facs[f]].phi;
    }
  return(TRUE);
}


#define MAX_LOG_N (1000)

mpfi_t log_n[MAX_LOG_N+1];

void simple_f_hat_even(unsigned int q, double eta, unsigned int n0, unsigned int n1, unsigned int N,double B,
		       unsigned int M, factor *factors, mpfi_c_t *chis, mpfi_ptr sqrt_q)
{
  debug;
  mpfi_c_t u_ieta,res,term;
  mpfi_c_init(u_ieta);mpfi_c_init(res);mpfi_c_init(term);
  double x;
  bool real_p;
  int prim=-1,prim1;
  unsigned int i,n,index,index1;
  mpfi_set(u_ieta->im,mpfi_pi_by_4);
  mpfi_mul_d(u_ieta->im,u_ieta->im,eta);
  for(index=0;index<factors[q].phi;index++)
    {
      if(!primitive(q,index,factors))
	continue;
      prim++;
      prim1=conj_j(prim,factors[q].phi,q);
      if(prim1<prim)
	continue;
      real_p=(prim==prim1);
      make_chis(q,index,factors,chis);
      for(i=1;i<q;i++)
	if(co_primes[i])
	  {
	    mpfi_set_ui(res->re,i);
	    mpfi_sqrt(res->re,res->re);
	    mpfi_c_div_i(chis[i],chis[i],res->re);
	  }
      mpfi_c_set(res,c_zero);
      for(x=n0*one_over_A,i=n0;i<n1;x+=one_over_A,i++)
	{
	  for(n=1;n<=M;n++)
	    {
	      if(!co_primes[n])
		continue;
	      if(n<=MAX_LOG_N)
		mpfi_set(u_ieta->re,log_n[n]);
	      else
		{
		  mpfi_set_ui(u_ieta->re,n);
		  mpfi_div(u_ieta->re,u_ieta->re,sqrt_q);
		  mpfi_log(u_ieta->re,u_ieta->re);
		}
	      mpfi_add_d(u_ieta->re,u_ieta->re,x);
	      G_even(term,u_ieta);
	      mpfi_c_mul(term,term,chis[n%q]);
	      mpfi_c_add(res,res,term);
	    }
	  printf("index %d G(%d)=",index,i);mpfi_c_print(res);
	}
      if(!real_p)
	{
	  for(i=1;i<q;i++)
	    if(co_primes[i])
	      mpfi_neg(chis[i]->im,chis[i]->im); // dirty conj
	  mpfi_c_set(res,c_zero);
	  for(x=n0*one_over_A,i=n0;i<n1;x+=one_over_A,i++)
	    {
	      for(n=1;n<=M;n++)
		{
		  if(!co_primes[n])
		    continue;
		  if(n<=MAX_LOG_N)
		    mpfi_set(u_ieta->re,log_n[n]);
		  else
		    {
		      mpfi_set_ui(u_ieta->re,n);
		      mpfi_div(u_ieta->re,u_ieta->re,sqrt_q);
		      mpfi_log(u_ieta->re,u_ieta->re);
		    }
		  mpfi_add_d(u_ieta->re,u_ieta->re,x);
		  G_even(term,u_ieta);
		  mpfi_c_mul(term,term,chis[n%q]);
		  mpfi_c_add(res,res,term);
		}
	      printf("index %d G(%d)=",index1,i);mpfi_c_print(res);
	    }
	}
    }
  mpfi_c_clear(u_ieta);mpfi_c_clear(term);mpfi_c_clear(res);
}

void FFT_f_hat_even(unsigned int q, double eta, unsigned int n0, unsigned int n1, unsigned int N,double B,unsigned int M,
		    factor *factors, mpfi_c_t *chis, mpfi_ptr sqrt_q)
{
  debug;
  unsigned int j,k,index;
  calc_polyf(TAYLOR_TERMS-2); // up to G^(TAYLOR_TERMS-1)(u) 
  for(index=0;index<factors[q].phi;index++)
    {
      make_chis(q,index,factors,chis);
      for(j=1;j<q;j++)
	if(co_primes[j])
	  mpfi_c_print(chis[j]);
      test_chis(q,factors,chis); // warning O(q^2)
    }
  /*
  for(k=0;k<TAYLOR_TERMS-3;k++)
    {
      for(j=0;j<k+2;j++)
	{
	  mpfi_out_str(stdout,10,0,polyf[k][j]);
	  printf("\n");
	}
      printf("\n");
    }
  */
}

int main(int argc, char **argv)
{
  int q1,n01,n11,prec,N1;
  unsigned int q,n0,n1,M,n,i,j,k,N;
  double eta,log2byq,cosetaetc,B;
  //mpfi_t x,delta,logdelta,pi_sindelta_by_q,err,log2,logq;
  //mpfi_t two_pi_by_B,eta_pi_by_4,mlog2byq,logpi,pibyq,pi;
  //mpfi_c_t res,term,outer_exp,exp_2_u_x,u_x;
  FILE *outfile,*facs_file;
  factor *factors;
  mpfi_c_t *chis;
  mpfi_t sqrt_q;


  if(argc!=9)
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
  for(i=1;i<q;i++)
    co_primes[i]=co_prime(i,q);

  eta=atof(argv[3]);
  if((eta<=-1.0)||(eta>=1.0))
    print_usage();

  n01=atoi(argv[4]);
  if(n01<0)
    print_usage();
  n0=(unsigned int) n01;

  n11=atoi(argv[5]);
  if(n11<=n01)
    print_usage();
  n1=(unsigned int) n11;

  N1=atoi(argv[6]);
  if((N1<(n1-1)*2)||((N1&63)!=0))
    print_usage();
  N=N1;
  B=N*one_over_A;
  
  facs_file=fopen(argv[7],"r");
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
  
  outfile=fopen(argv[8],"wb");
  if(!outfile)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[5]);
      exit(1);
    }

  fwrite(&q,sizeof(unsigned int),1,outfile);
  fwrite(&eta,sizeof(double),1,outfile);
  fwrite(&n0,sizeof(unsigned int),1,outfile);
  fwrite(&n1,sizeof(unsigned int),1,outfile);
  fwrite(&N,sizeof(unsigned int),1,outfile);
  if(!(chis=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*q)))
    {
      printf("Failed to allocate memory for chis.Exiting\n");
      exit(1);
    }
  for(i=0;i<q;i++)
    mpfi_c_init(chis[i]);
  mpfi_init(theta);
  mpfi_init(theta1);
  mpfi_c_init(chi);
  mpfi_c_init(chi_neg);
  mpfi_c_init(mpfi_c_small);
  mpfi_c_init(G_even_tmp);
  mpfi_set_d(mpfi_c_small->re,LN_SMALL);
  mpfi_exp(mpfi_c_small->re,mpfi_c_small->re);
  mpfi_neg(mpfi_c_small->im,mpfi_c_small->re);
  mpfi_put(mpfi_c_small->re,mpfi_c_small->im);
  mpfi_set(mpfi_c_small->im,mpfi_c_small->re);
  printf("small=");mpfi_c_print(mpfi_c_small);
  mpfi_init(sqrt_q);
  mpfi_set_ui(sqrt_q,q);
  mpfi_sqrt(sqrt_q,sqrt_q);  

  for(n=1;n<=MAX_LOG_N;n++)
    {
      mpfi_init(log_n[n]);
      mpfi_set_ui(log_n[n],n);
      mpfi_div(log_n[n],log_n[n],sqrt_q);
      mpfi_log(log_n[n],log_n[n]);
    }

  M=calc_M(2.0*n*M_PI/B,log((double) 2.0)-log((double) q)*0.25,cos(M_PI*eta/2.0)*M_PI/q);
  printf("We are using %d terms in f_hat sum.\n",M);
  
  if(M<MIN_M_FOR_FFT)
    simple_f_hat_even(q,eta,n0,n1,N,B,M,factors,chis,sqrt_q);
  else
    FFT_f_hat_even(q,eta,n0,n1,N,B,M,factors,chis,sqrt_q);



  exit(0);
  /*
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
  write_mpfi(outfile,two_pi_by_B);
  mpfi_init(pibyq);
  mpfi_div_ui(pibyq,pi,q);
  mpfi_init(logpi);
  mpfi_log(logpi,pi);
  mpfi_mul_d(delta,pi,(1.0-fabs(eta))/2.0);
  write_mpfi(outfile,delta);
  mpfi_init(logdelta);
  mpfi_log(logdelta,delta);
  write_mpfi(outfile,logdelta);
  mpfi_mul_d(eta_pi_by_4,mpfi_pi_by_4,eta);
  write_mpfi(outfile,eta_pi_by_4);
  mpfi_init(pi_sindelta_by_q);
  mpfi_sin(pi_sindelta_by_q,delta);
  mpfi_mul(pi_sindelta_by_q,pi_sindelta_by_q,pi);
  mpfi_div_ui(pi_sindelta_by_q,pi_sindelta_by_q,q);
  write_mpfi(outfile,pi_sindelta_by_q);
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
  mpfi_c_init(res);
  mpfi_c_init(outer_exp);
  mpfi_c_init(inner_exp);
  mpfi_c_init(u_x);
  mpfi_c_init(exp_2_u_x);
  mpfi_init(log_n[0]);
  mpfi_init(log_n[1]);
  mpfi_set_ui(log_n[1],0);
  for(n=1;n<=MAX_LOG_N;n++)
    {
      mpfi_init(log_n[n]);
      mpfi_set_ui(log_n[n],n);
      mpfi_log(log_n[n],log_n[n]);
    }
  mpfi_init(log2);
  mpfi_init(logq);
  mpfi_set_d(log2,2.0);
  mpfi_log(log2,log2);
  write_mpfi(outfile,log2);
  mpfi_set_ui(logq,q);
  mpfi_log(logq,logq);
  write_mpfi(outfile,logq);
  mpfi_mul_d(logq,logq,0.75);
  mpfi_init(mlog2byq);
  mpfi_sub(mlog2byq,log2,logq);
  mpfi_set(u_x->im,eta_pi_by_4);
  for(n=n0;n<n1;n++)
    {
      mpfi_mul_ui(x,two_pi_by_B,n);
      mpfi_set(u_x->re,x);
      mpfi_c_mul_d(exp_2_u_x,u_x,2.0);
      mpfi_c_exp(exp_2_u_x,exp_2_u_x);
      mpfi_c_mul_i(exp_2_u_x,exp_2_u_x,pibyq);
      //printf("exp(2U(x))*pi/q=");mpfi_c_print(exp_2_u_x);
      mpfi_c_mul_d(outer_exp,u_x,1.5);
      mpfi_add(outer_exp->re,outer_exp->re,log2);
      mpfi_sub(outer_exp->re,outer_exp->re,logq);
      //printf("M=%d\n",M);
      fwrite(&M,sizeof(unsigned int),1,outfile);
      F_hat_err(err,x,pi_sindelta_by_q,mlog2byq,M);
      write_mpfi(outfile,err);
      //printf("err=");mpfi_print(err);

      for(i=1;i<min(q,M+1);i++)
	if(co_prime(i,q))
	  {
	    mpfi_set_ui(res->re,0);
	    mpfi_set_ui(res->im,0);
	    for(k=i;k<=M;k+=q)
	      {
		F_hat_o_term(term,outer_exp,exp_2_u_x,k);
		mpfi_c_inc(res,term);
	      }
	    //printf("FFT=");mpfi_c_print(res);
	    write_mpfi_c(outfile,res);
	  }
    }
  fclose(outfile);
  */
}



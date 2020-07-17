/*

File: G1.6.c

Created: 15 March 2011

Version: 1.6

Last Modified: 

Dialect: C

Requires: GMP
          MPFR
          MPFI

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

V 1.0 Initial implementation
V 1.1 Control over maximum size of h added
      Calculation of G(1) added
      mpfi_c routines moved to mpfi_c.c
V 1.2 Revised Taylor Expansion used and taylor error added
V 1.3 Now takes zeros file in standard form
v 1.4 Taylor Error Fixed. Lambda now entered as 2^(-n)
v 1.5 rewritten for uber-accurate zeros from windowed zeta
      now uses mpfi_c.h
v 1.6 Improved Taylor error


Build instructions: gcc -oG1.5 G1.5.c -O2 -lmpfi

By: DJ Platt
    Bristol University

Copyright 2008,2009,2010,2011.

The author is funded by the UK ESPRC. */

#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "time.h"
#include "../includes/mpfi_c.h"
#include "../includes/pi_x.h"
#include "../windowed/win_zeta.h"

double H_MAX;       /* absolute upper limit for h_max otherwise taylor won't converge */
#define SUCCESS (0)



/* global variables for commonly used values initialised in main */

mpfr_t lambda;
mpfi_t lam_sqr,lam_sqr_by_2,lnx;
//mpz_t x;
long int taylor_n;
mpfi_t h_max;
mpfi_c_t F_tmp1,F_tmp2;

inline void next_rho(mpfi_t del_t, FILE *infile)
{
  in_bytes(del_t,infile);
}

/* F(s)=exp(lambda^2*s^2/2)*x^s/s */
void F (mpfi_c_ptr res,mpfi_c_ptr s)
{
    mpfi_c_mul(F_tmp1,s,s);
    mpfi_c_mul_i(F_tmp1,F_tmp1,lam_sqr_by_2);
    mpfi_c_mul_i(F_tmp2,s,lnx);
    mpfi_c_add(F_tmp1,F_tmp1,F_tmp2);
    mpfi_c_exp(F_tmp1,F_tmp1);
    mpfi_c_div(res,F_tmp1,s);
};


/* global variables used by G
   initialised in init_vectors()
   destroyed in free_vectors()
*/

mpfi_c_t F_s,zero,*ihk_n,k,k_n,one,*ihk_n_sum_pos,*ihk_n_sum_neg;
mpfi_c_t *s_minus_n,mit,msum,nsum_pos,nsum_neg,nsum,taylor_err;
mpfi_t *lam_sqr_by_2_n,Ea,Eb,R,one_R;
mpz_t *nbang;

/* adds the real taylor error at s and h,-h to err */
void calc_taylor_error(mpfi_c_ptr ans, double t0, double H, unsigned long int taylor_n)
{
  mpfi_mul_d(Ea,lam_sqr_by_2,H);
  mpfi_mul_d(Ea,Ea,H);
  mpfi_log(Ea,Ea);
  mpfi_mul_ui(Ea,Ea,taylor_n/2);
  mpfi_exp(Ea,Ea);
  mpfi_div_z(Ea,Ea,nbang[taylor_n/2]); // =(lam^2*H^2/2)^(n/2)/(n/2)!
  //mpfi_print_str("Ea=",Ea);
  mpfi_set_d(R,H);
  mpfi_div_d(R,R,t0); // R passes for H/|s|
  mpfi_set_ui(one_R,1);
  mpfi_sub(one_R,one_R,R); // 1-R
  mpfi_log(R,R);
  mpfi_mul_ui(R,R,taylor_n);
  mpfi_exp(R,R);
  mpfi_div(Eb,R,one_R);    // R^n/(1-R)
  //mpfi_print_str("Eb=",Eb);
  mpfi_div(taylor_err->re,Ea,one_R);
  mpfi_add(taylor_err->re,taylor_err->re,Eb);
  mpfi_mul(taylor_err->im,Ea,Eb);
  mpfi_add(taylor_err->re,taylor_err->re,taylor_err->im); // EaEb+Eb+Ea/(1-R)
  mpfi_div_2ui(Ea,lam_sqr_by_2,2); // lam^2/8
  mpfi_mul_d(Eb,lam_sqr_by_2,t0);
  mpfi_mul_d(Eb,lam_sqr_by_2,t0);  // lam^2*t^2/2
  mpfi_sub(Ea,Ea,Eb);
  mpfi_div_2ui(Eb,lnx,1);          // log(x)/2
  mpfi_add(Ea,Ea,Eb);
  mpfi_exp(Ea,Ea);                 // sqrt(x)*exp(lam^2/8-lam^2*t^2/2)
  mpfi_div_d(Ea,Ea,t0);            // |F(s0)|
  mpfi_mul_d(Ea,Ea,H*2);           
  mpfi_mul(taylor_err->re,taylor_err->re,Ea); //2H|F(s0)|(EaEb+Eb+Ea/(1-R))
  mpfi_neg(taylor_err->im,taylor_err->re);
  mpfi_put(taylor_err->re,taylor_err->im);
  mpfi_set(taylor_err->im,taylor_err->re);
  //mpfi_c_print_str("Taylor error now=",taylor_err);
};



int G(mpfi_c_ptr ans, mpfi_c_ptr s, mpfi_ptr h)
{

/* returns the taylor estimate of G(s+ih)-G(s-ih) using
   taylor_n iterations.
*/
  	long unsigned int n,m;

	//mpfi_c_print_str("In G with s=",s);
	//mpfi_print_str("and h=",h);
	F(F_s,s);	// calculate F(s)
	//mpfi_c_print_str("F(s)=",F_s);
	

// set ihk_n[1]=i*h*(lambda^2*s+ln(x))
// set ihk_n_sum_neg[n]=sum 0..n ihk_n[n]
// set ihk_n_sum_pos[n]=sum 0..n (-1)^n*ihk_n[n]
	mpfi_c_set(k,s);
	mpfi_c_mul_i(k,k,lam_sqr);
	mpfi_c_add_re(k,k,lnx);
	//mpfi_c_print_str("ihk_n[1]=",ihk_n[1]);
	mpfi_c_mul_i(ihk_n[1],k,h);
	//mpfi_c_print_str("ihk_n[1]=",ihk_n[1]);
	mpfi_c_muli(ihk_n[1]);
	//mpfi_c_print_str("ihk_n[1]=",ihk_n[1]);

// set ihk_n[2..taylor_n*2]=(ihk)^n/n!

	mpfi_c_add(ihk_n_sum_neg[1],ihk_n[1],ihk_n_sum_neg[0]);
	mpfi_c_sub(ihk_n_sum_pos[1],ihk_n_sum_pos[0],ihk_n[1]);
	
	for(n=2;n<=taylor_n+taylor_n;n++)
	{
	    mpfi_c_mul(ihk_n[n],ihk_n[n-1],ihk_n[1]);
	    mpfi_c_div_ui(ihk_n[n],ihk_n[n],n);
	    mpfi_c_add(ihk_n_sum_neg[n],ihk_n_sum_neg[n-1],ihk_n[n]);
	    if(n%2==0)
		mpfi_c_add(ihk_n_sum_pos[n],ihk_n_sum_pos[n-1],ihk_n[n]);
		else
		mpfi_c_sub(ihk_n_sum_pos[n],ihk_n_sum_pos[n-1],ihk_n[n]);
	};

	

// set k_n=k^1=k;

	mpfi_c_set(k_n,k);
	
	for(n=1;n<=taylor_n+taylor_n;n++)
	    mpfi_c_div(s_minus_n[n],s_minus_n[n-1],s);  // set up (1/s)^n

	mpfi_c_set(nsum_pos,zero);
	mpfi_c_set(nsum_neg,zero);
	
	for(n=0;n<=taylor_n;n++)
	{
	    mpfi_c_set(msum,zero);	    
	    for(m=0;m<=(n>>1);m++)
	    {
		mpfi_c_mul_i(mit,s_minus_n[n-m-m],lam_sqr_by_2_n[m]);
		mpfi_c_div_z(mit,mit,nbang[m]);
		mpfi_c_add(msum,msum,mit);		
	    };
	    mpfi_c_mul(nsum,msum,ihk_n_sum_pos[n]);
	    mpfi_c_mul_z(nsum,nsum,nbang[n]);
	    mpfi_c_div(nsum,nsum,k_n);
	    mpfi_c_add(nsum_pos,nsum_pos,nsum);
	    mpfi_c_mul(nsum,msum,ihk_n_sum_neg[n]);
	    mpfi_c_mul_z(nsum,nsum,nbang[n]);
	    mpfi_c_div(nsum,nsum,k_n);
	    mpfi_c_add(nsum_neg,nsum_neg,nsum);
	    


	    mpfi_c_mul(k_n,k_n,k);
	};

	for(n=taylor_n+1;n<=taylor_n+taylor_n;n++)
	{
	    mpfi_c_set(msum,zero);
	    for(m=(n+1)>>1;m<=taylor_n;m++)
	    {
		mpfi_c_mul_i(mit,s_minus_n[taylor_n+n-m-m],lam_sqr_by_2_n[m-(taylor_n>>1)]);
		mpfi_c_div_z(mit,mit,nbang[m-(taylor_n>>1)]);
		mpfi_c_add(msum,msum,mit);				
	    };

	    
	    mpfi_c_mul_z(msum,msum,nbang[n]);
	    mpfi_c_div(msum,msum,k_n);
	    mpfi_c_mul(nsum,msum,ihk_n_sum_pos[n]);
	    mpfi_c_add(nsum_pos,nsum_pos,nsum);

	    mpfi_c_mul(nsum,msum,ihk_n_sum_neg[n]);
	 
	    mpfi_c_add(nsum_neg,nsum_neg,nsum);
	    


	    mpfi_c_mul(k_n,k_n,k);
	};
	
	//mpfi_c_print_str("msum=",msum);
	//mpfi_c_print_str("ihk[1]=",ihk_n[1]);
	mpfi_c_neg(msum,ihk_n[1]);       // msum = -ihk
	
	//mpfi_c_print_str("msum=",msum);
	mpfi_c_exp(msum,msum);		 // msum = exp(-ihk)
	
	mpfi_c_mul(msum,msum,nsum_neg);  // msum = exp(-ihk)*sum
	
	mpfi_c_exp(ans,ihk_n[1]);        // ans  = exp(ihk)
	
	mpfi_c_mul(ans,ans,nsum_pos);    // ans  = exp(ihk)*sum
	
	mpfi_c_sub(ans,msum,ans);
	
	mpfi_c_mul(ans,ans,F_s);

	mpfi_c_neg(ans,ans);
	
	mpfi_c_add(ans,ans,taylor_err);
	//mpfi_c_print_str("G returning ",ans);

}; /* G */

int init_vectors()
{
    
    long int n;

    mpfi_c_init(taylor_err);
    mpfi_init(R);
    mpfi_init(one_R);
    mpfi_init(Ea);
    mpfi_init(Eb);
    mpfi_c_init(F_s);
    mpfi_c_init(zero);
    mpfi_c_init(k);
    mpfi_c_init(k_n);
    mpfi_c_init(one);
    mpfi_c_init(mit);
    mpfi_c_init(msum);
    mpfi_c_init(nsum_pos);
    mpfi_c_init(nsum_neg);
    mpfi_c_init(nsum);
    
    mpfi_set_ui(one->re,1);
    mpfi_set_ui(one->im,0);  // one <- 1+0i

    
    mpfi_set_ui(zero->re,0);
    mpfi_set_ui(zero->im,0);


    nbang=malloc((taylor_n+taylor_n+1) * sizeof(mpz_t));
    if(nbang==NULL)
    {
	printf("Error allocating memory for nbang.\n");
	return(FAILURE);
    };    
    

    mpz_init(nbang[0]);
    mpz_set_ui(nbang[0],1);            /* nbang[0]<-0! */
    for(n=1;n<=taylor_n+taylor_n;n++)
    {
	mpz_init(nbang[n]);
	mpz_mul_ui(nbang[n],nbang[n-1],n); /* nbang[n]<-n! */
    };


    ihk_n=malloc((taylor_n+taylor_n+1) * sizeof(mpfi_c_t));
    if(ihk_n==NULL)
    {
	printf("Error allocating memory for ihk_n.\n");
	return(FAILURE);
    };
    for(n=0;n<=taylor_n+taylor_n;n++)
	mpfi_c_init(ihk_n[n]);
    mpfi_c_set(ihk_n[0],one);      // ihk_n[0] <- 1+0i=(ihk)^0


    ihk_n_sum_pos=malloc((taylor_n+taylor_n+1) * sizeof(mpfi_c_t));
    if(ihk_n_sum_pos==NULL)
    {
	printf("Error allocating memory for ihk_n_sum_pos.\n");
	return(FAILURE);
    };
    for(n=0;n<=taylor_n+taylor_n;n++)
	mpfi_c_init(ihk_n_sum_pos[n]);
    mpfi_c_set(ihk_n_sum_pos[0],one);

    ihk_n_sum_neg=malloc((taylor_n+taylor_n+1) * sizeof(mpfi_c_t));
    if(ihk_n_sum_neg==NULL)
    {
	printf("Error allocating memory for ihk_n_sum_neg.\n");
	return(FAILURE);
    };
    for(n=0;n<=taylor_n+taylor_n;n++)
	mpfi_c_init(ihk_n_sum_neg[n]);
    mpfi_c_set(ihk_n_sum_neg[0],one);

    lam_sqr_by_2_n=malloc((taylor_n+1) * sizeof(mpfi_t));
    if(lam_sqr_by_2_n==NULL)
    {
	printf("Error allocating memory for lam_sqr_by_2_n.\n");
	return(FAILURE);
    };
    for(n=0;n<=taylor_n;n++)
	mpfi_init(lam_sqr_by_2_n[n]);
    mpfi_set_ui(lam_sqr_by_2_n[0],1);
    mpfi_set(lam_sqr_by_2_n[1],lam_sqr_by_2);
    for(n=2;n<=taylor_n;n++)
	mpfi_mul(lam_sqr_by_2_n[n],lam_sqr_by_2_n[n-1],lam_sqr_by_2);

    s_minus_n=malloc((taylor_n+taylor_n+1) * sizeof(mpfi_c_t));
    if(s_minus_n==NULL)
    {
	printf("Error allocating memory for s_minus_n.\n");
	return(FAILURE);
    };
    for(n=0;n<=taylor_n+taylor_n;n++)
	mpfi_c_init(s_minus_n[n]);
    mpfi_c_set(s_minus_n[0],one);



    return(SUCCESS);
};



int main(int argc, char **argv)
{
  FILE *infile;
  if(argc!=5)
    {
      printf("Usage G1.5 <prec> <taylor_n> <h max> <infile>\n");
      printf("h_max must be in (0,0.5]\n");
      printf("taylor_n must be >0 and even.\n");
      exit(0);
    }
  if(!(infile=fopen(argv[4],"rb")))
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[4]);
      exit(0);
    }
  mpfi_c_setup(atoi(argv[1]));
  H_MAX=atof(argv[3]);
  if((H_MAX<=0.0)||(H_MAX>0.5))
    {
      printf("Usage G1.5 <prec> <taylor_n> <h max> <infile>\n");
      printf("h_max must be in (0,0.5]\n");
      printf("taylor_n must be >0 and even.\n");
      exit(0);
    }
  mpfi_c_init(F_tmp1);
  mpfi_c_init(F_tmp2);
  mpfi_init(lam_sqr);
  mpfi_init(lam_sqr_by_2);
  mpfi_set_d(lam_sqr,LAMBDA);
  mpfi_sqr(lam_sqr,lam_sqr);
  mpfi_init(lnx);
  mpfi_set_128(lnx,calc_x());
  printf("x=");print_bigint(calc_x());
  mpfi_log(lnx,lnx);
  mpfi_init(h_max);
  mpfi_set_d(h_max,H_MAX);

  mpfi_div_2ui(lam_sqr_by_2,lam_sqr,1);
  taylor_n=atoi(argv[2]);
  if((taylor_n<2)||((taylor_n&1)==1))
    {
      printf("Usage G1.5 <prec> <taylor_n> <h max> <infile>\n");
      printf("h_max must be in (0,0.5]\n");
      printf("taylor_n must be >0 and even.\n");
      exit(0);
    }
  init_vectors();
  printf(" Processing file %s ",argv[4]);
  printf("Running with lambda=2^%f H_MAX=%f taylor_n=%ld\n",log(LAMBDA)/log(2.0),H_MAX,taylor_n);
  mpfi_c_t res,res1,res2,res3,s;
  mpfi_t h,t,del_t,del_h;
  long int num_its,it,z;
  double st[2],t0,t1;
  long int zs[2],m;
  mpfi_c_init(res);
  mpfi_c_init(res1);
  mpfi_c_init(res2);
  mpfi_c_init(res3);
  mpfi_c_zero(res3);
  mpfi_c_init(s);
  mpfi_init(h);
  mpfi_init(t);
  mpfi_init(del_t);
  mpfi_init(del_h);
  mpfi_c_zero(res);
  mpfi_set_d(s->re,0.5);
  init_in_bytes();
  fread(&num_its,sizeof(long int),1,infile);
  //printf("Doing %ld iterations.\n",num_its);
  int calced=(1==0);
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile);
      fread(&zs[0],sizeof(long int),1,infile);
      if(st[0]==0.0)
	continue;
      if(!calced)
	{
	  calced=(1==1);
	  t0=st[0];
	  //n=zs[0]+1;
	  calc_taylor_error(taylor_err,t0,H_MAX,taylor_n);
	}
      fread(&zs[1],sizeof(long int),1,infile);
      //printf("Processing zero %ld to %ld=%ld in total.\nt0=%f\n",zs[0]+1,zs[1],zs[1]-zs[0],t0);mpfi_c_print_str("Taylor error=",taylor_err);exit(0);
      mpfi_set_d(s->im,st[0]);
      mpfi_set(t,s->im);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  next_rho(del_t,infile);
          if(mpfi_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  mpfi_add(t,t,del_t); // exact, used to avoid accumulation of +/- 2^-(OP_ACC+1)
	  mpfi_add(del_t,del_t,pm1); // not exact
	  //mpfi_print_str("t   =",t);
	  //mpfi_print_str("del_t=",del_t);
	  if(mpfi_cmp_d(del_t,H_MAX*2.0)>0) // can't get there in one step
	    {
	      mpfi_add_d(s->im,s->im,H_MAX);
	      G(res1,s,h_max);
	      mpfi_sub_d(del_t,del_t,H_MAX*2.0);
	      while(mpfi_cmp_d(del_t,H_MAX*2.0)>0) // do big steps
		{
		  mpfi_add_d(s->im,s->im,H_MAX*2.0);
		  G(res2,s,h_max);
		  mpfi_c_add(res1,res1,res2);
		  mpfi_sub_d(del_t,del_t,H_MAX*2.0);
		}
	      mpfi_add_d(s->im,s->im,H_MAX); // now do a baby step
	      //mpfi_c_print_str("Reached to within 2*H_MAX with s=",s);
	      //mpfi_print_str("and with del_t=",del_t);
	      mpfi_div_2ui(del_h,del_t,1);
	      mpfi_add(s->im,s->im,del_h);
	      G(res2,s,del_h);
	      mpfi_c_add(res1,res1,res2);
	      mpfi_set(s->im,t);
	    }
	  else // can get there in one baby step
	    {
	      mpfi_div_2ui(del_h,del_t,1);
	      mpfi_add(s->im,s->im,del_h);
	      G(res1,s,del_h);
	      mpfi_set(s->im,t);
	    }

	  //printf("zero number %ld gave",z);
	  //mpfi_c_print_str("",res1);

	  mpfi_c_add(res3,res3,res1);
	  mpfi_c_mul_ui(res1,res1,z);
	  mpfi_c_add(res,res,res1);
	  //exit(0);
	}
      // now do the G(end)-G(s_last)
      //printf("st[1]=%12.0f\n",st[1]);
      //mpfi_c_print_str("s=",s);
      mpfi_sub_d(del_t,s->im,st[1]);
      mpfi_neg(del_t,del_t);
      //mpfi_print_str("Final delta t = ",del_t);
      mpfi_set_d(t,st[1]);
      if(mpfi_cmp_d(del_t,H_MAX*2.0)>0)
	{
	  mpfi_add_d(s->im,s->im,H_MAX);
	  G(res1,s,h_max);
	  mpfi_sub_d(del_t,del_t,H_MAX*2.0);
	  while(mpfi_cmp_d(del_t,H_MAX*2.0)>0)
	    {
	      mpfi_add_d(s->im,s->im,H_MAX*2.0);
	      G(res2,s,h_max);
	      mpfi_c_add(res1,res1,res2);
	      mpfi_sub_d(del_t,del_t,H_MAX*2.0);
	    }
	  mpfi_add_d(s->im,s->im,H_MAX);
	  mpfi_div_2ui(del_h,del_t,1);
	  mpfi_add(s->im,s->im,del_h);
	  G(res2,s,del_h);
	  mpfi_c_add(res1,res1,res2);
	  mpfi_set(s->im,t);
	}
      else
	{
	  mpfi_div_2ui(del_h,del_t,1);
	  mpfi_add(s->im,s->im,del_h);
	  G(res1,s,del_h);
	  mpfi_set(s->im,t);
	}
      mpfi_c_add(res3,res3,res1);
      //mpfi_c_print_str("res3 now=",res3);
      mpfi_c_mul_ui(res1,res1,z);
      mpfi_c_add(res,res,res1);
      //mpfi_c_print_str("res now =",res);
      //exit(0);
    }
  t1=st[1];
  m=zs[1];
  //printf("(%ld)*G(%f)-sum_{k=%ld}^%ld G(rho_k)-%ld*G(%f)=",m+1,t1,n,m,n,t0);
  mpfi_c_print_str("",res);
  //printf("G(%f)-G(%f)=",t1,t0);
  mpfi_c_print_str("",res3);
  return(0);
}

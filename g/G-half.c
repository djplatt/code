/*

File: G-half.c

Created: 4 March 2008

Version: 1.0

Last Modified: 4 March 2008

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

Build instructions: gcc -oG1.1 G-half.c -O2 -lmpfi -lmpfr -lgmp

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "mpfi_c.c"

#define SUCCESS 0
#define QUIT_CODE 0
#define FAILURE 1
#define TRUE (1==1)
#define FALSE (0==1)
#define H_MAX 0.5       /* absolute upper limit for h_max */
#define max_height 100.0

/* global varibale for commonly used values */

mpfr_t lambda;
mpfi_t lam_sqr,lam_sqr_by_2,lnx;
mpz_t x;
int taylor_n;


/* temporary variables for mul/div/exp/F
   declared globally to avoid having to initialise
   each time mul etc. gets called */


/* the imaginary parts of the zeros will be kept in here */



void F (mpfi_c_ptr res,mpfi_c_ptr s)
{
//    printf("s is ");mpfi_c_print(s);
    mpfi_c_mul(F_tmp1,s,s);
//    printf("s^2 is ");mpfi_c_print(F_tmp1);
    mpfi_c_mul_i(F_tmp1,F_tmp1,lam_sqr_by_2);
//    printf("(s*lambda)^2/2 is ");mpfi_c_print(F_tmp1);
    mpfi_c_mul_i(F_tmp2,s,lnx);
//    printf("sln(x) is ");mpfi_c_print(F_tmp2);
    mpfi_c_add(F_tmp1,F_tmp1,F_tmp2);
//    printf("(s*lambda)^2/2+sln(x) is ");mpfi_c_print(F_tmp1);
    mpfi_c_exp(F_tmp1,F_tmp1);
//    printf("exp(...) is ");mpfi_c_print(F_tmp1);
    mpfi_c_div(res,F_tmp1,s);
};


mpfi_c_t ks,F_s,n_sum_pos;
mpfi_c_t n_sum_neg,ks_sqr,c_tmp1,zero;
mpz_t *nbang;
mpfi_c_t *ihk_n;
mpfi_c_t *ks_2n_kks;        // k^(2n+2)s^(2n+1)
mpfi_c_t *m_sum2;    // sigma (lambda^2*s^2/2)^n/n!
mpfi_c_t *m_sum1_pos;
mpfi_c_t *m_sum1_neg;
mpfr_t h_max;





/* this sets err to be the lower and upper bound for the real integral error,
   i.e. the error introduced by taking G(rho_max)=0
   not yet implemented so adds 0+0i (nop)
*/
void calc_integral_error(mpfi_ptr err)
{
    
};

/* this adds to err the lower and upper bounds for the real error introduced
   by stopping the sum over rho of G at rho_max
   not yet implemented so adds [0,0] (nop)
*/
void calc_trunc_error(mpfi_ptr err)
{
};

/* this adds the real taylor error at s,h
   and s,-h to err
   not yet implemented so adds [0,0] (nop)
*/

/* adds the real taylor error at s and h,-h to err */
void calc_taylor_error(mpfi_c_ptr err,mpfi_c_ptr s,mpfi_ptr h)
{
};





int G(mpfi_c_ptr ans, mpfi_c_ptr s, mpfi_ptr h)
{

/* returns the taylor estimate of G(s+ih)-G(s-ih) using
   taylor_n iterations.
*/
  	int n;

	mpfi_c_mul_i(ks,s,lam_sqr);
	mpfi_c_add_re(ks,ks,lnx);   // ks <- k
	mpfi_c_set(ks_2n_kks[0],ks); // ks_2n_kks <- k
	mpfi_c_mul_i(ihk_n[1],ks,h);     // ihk_n[1] <- hk
	mpfi_c_mul(ks,ks,s);        // ks <- ks
	mpfi_c_muli(ihk_n[1]);           // ihk_n[1] <- ihk

        for(n=2;n<taylor_n+taylor_n+2;n++)
	    mpfi_c_mul(ihk_n[n],ihk_n[n-1],ihk_n[1]);  // ihk[n]<-(ihk)^n
//	for(n=0;n<taylor_n+taylor_n+2;n++)
//            {printf("ihk^%d is ",n);mpfi_c_print(ihk_n[n]);}
	for(n=1;n<=taylor_n;n++)
	{
	    mpfi_c_div_z(ans,ihk_n[n+n-1],nbang[n+n-1]);
	    mpfi_c_div_z(c_tmp1,ihk_n[n+n],nbang[n+n]);
	    mpfi_c_add(m_sum1_pos[n],m_sum1_pos[n-1],c_tmp1);
	    mpfi_c_add(m_sum1_neg[n],m_sum1_neg[n-1],c_tmp1);
	    mpfi_c_sub(m_sum1_pos[n],m_sum1_pos[n],ans);
	    mpfi_c_add(m_sum1_neg[n],m_sum1_neg[n],ans);
	};

	mpfi_c_mul(ks_sqr,ks,ks);
	mpfi_c_mul(ks_2n_kks[0],ks_2n_kks[0],ks);   // ks_2n_kks[0] <- kks
	for(n=1;n<=taylor_n;n++)
	    mpfi_c_mul(ks_2n_kks[n],ks_2n_kks[n-1],ks_sqr); // ks_2n_kks[n]
                                                        // <- (ks)^(2n)*kks
//	for(n=0;n<=taylor_n;n++)
//	    {printf("k^%ds%d is ",2*n+2,2*n+1);mpfi_c_print(ks_2n_kks[n]);}

	mpfi_c_mul(m_sum2[1],s,s);                      // m_sum2[1]<-s^2
	mpfi_c_mul_i(m_sum2[1],m_sum2[1],lam_sqr_by_2); // m_sum2[1]<-l^2s^2/2
	for(n=2;n<=taylor_n;n++)
	    mpfi_c_mul(m_sum2[n],m_sum2[n-1],m_sum2[1]);  // m_sum2[n] 
	                                               // <- (l^2s^2/2)^n
	for(n=2;n<=taylor_n;n++)
	    mpfi_c_div_z(m_sum2[n],m_sum2[n],nbang[n]);
	    
	for(n=1;n<=taylor_n;n++)
	    mpfi_c_add(m_sum2[n],m_sum2[n],m_sum2[n-1]);

	mpfi_c_set(n_sum_pos,zero);
	mpfi_c_set(n_sum_neg,zero);
	for(n=0;n<=taylor_n;n++)
	{
	    
	    mpfi_c_add_ui(ans,ks,n+n+1);            // ans<-ks+2n+1
	    mpfi_c_mul_z(ans,ans,nbang[n+n]);     // ans<-(2n)!(ks+2n+1)
//	    printf("(2n)!(ks+2n+1) is ");mpfi_c_print(ans);
	    mpfi_c_mul(c_tmp1,ans,m_sum1_pos[n]);// c_tmp1<-(2n)!(ks+2n+1)M1+
	    mpfi_c_mul(ans,ans,m_sum1_neg[n]);  // ans<-(2n)!(ks+2n+1)M1-
	    mpfi_c_sub(c_tmp1,c_tmp1,ihk_n[n+n+1]); // c_tmp1<-~-(ihk)^(2n+1)
	    mpfi_c_add(ans,ans,ihk_n[n+n+1]);   // ans<-~+(ihk)^(2n+1)
	    mpfi_c_mul(c_tmp1,c_tmp1,m_sum2[n]);    
	    mpfi_c_mul(ans,ans,m_sum2[n]);
	    mpfi_c_div(c_tmp1,c_tmp1,ks_2n_kks[n]);
	    mpfi_c_div(ans,ans,ks_2n_kks[n]);
	    mpfi_c_add(n_sum_pos,n_sum_pos,c_tmp1);
	    mpfi_c_add(n_sum_neg,n_sum_neg,ans);
//	    printf("iteration with n=%d\n",n);
//	    printf("m_sum1_pos is ");mpfi_c_print(m_sum1_pos[n]);
//	    printf("m_sum1_neg is ");mpfi_c_print(m_sum1_neg[n]);
//	    printf("m_sum2 is ");mpfi_c_print(m_sum2[n]);
//	    printf("pos adds ");mpfi_c_print(c_tmp1);
//	    printf("neg adds ");mpfi_c_print(ans);
//	    printf("n_sum_pos is ");mpfi_c_print(n_sum_pos);
//	    printf("n_sum_neg is ");mpfi_c_print(n_sum_neg);
	};

	mpfi_c_exp(ans,ihk_n[1]);
//	printf("exp(ihk) is ");mpfi_c_print(ans);
	mpfi_c_mul(ans,ans,n_sum_pos);
//	printf("e(ihk)*n_sum_pos is ");mpfi_c_print(ans);
	mpfi_c_neg(c_tmp1,ihk_n[1]);
	mpfi_c_exp(c_tmp1,c_tmp1);
//	printf("exp(-ihk) is ");mpfi_c_print(c_tmp1);
	mpfi_c_mul(c_tmp1,c_tmp1,n_sum_neg);
//	printf("e(-ihk)*n_sum_neg is ");mpfi_c_print(c_tmp1);
	mpfi_c_sub(ans,ans,c_tmp1);
	F(F_s,s);
	mpfi_c_mul(ans,F_s,ans);

//	printf("s was ");mpfi_c_print(s);
//	printf("h was ");mpfi_print(h);
//	printf("result ");mpfi_c_print(ans);

}; /* G */



int g_half()
{
    unsigned long long int i;
    unsigned long int n,m;
    mpfi_t half,h,tmp,h_max_i;
    mpfr_t h_fr,two_h_max;
    mpfi_c_t s,g_rho_delta,one;
    mpfi_c_t g_rho_delta_sum;


    mpfi_init(half);


/* assign 1/2 to half */
    if(mpfi_set_d(half,0.5)!=MPFI_FLAGS_BOTH_ENDPOINTS_EXACT)
    {
	printf("Inconceivable error setting up 1/2.\n");
	return(FAILURE);
    };

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
    ihk_n=malloc((taylor_n+taylor_n+2) * sizeof(mpfi_c_t));
    if(ihk_n==NULL)
    {
	printf("Error allocating memory for ihk_n.\n");
	return(FAILURE);
    };
//
    for(n=0;n<taylor_n+taylor_n+2;n++)
	mpfi_c_init(ihk_n[n]);
    mpfi_c_init(one);
    mpfi_set_ui(one->re,1);
    mpfi_set_ui(one->im,0);
    mpfi_c_set(ihk_n[0],one);      // ihk_n[0] <- 1+0i=(ihk)^0

    ks_2n_kks=malloc((taylor_n+1) * sizeof(mpfi_c_t));
    if(ks_2n_kks==NULL)
    {
	printf("Error allocating memory for ks_2n_kks.\n");
	return(FAILURE);
    };
    for(n=0;n<=taylor_n;n++)
	mpfi_c_init(ks_2n_kks[n]);

    m_sum2=malloc((taylor_n+1) * sizeof(mpfi_c_t));
    if(m_sum2==NULL)
    {
	printf("Error allocating memory for ls2_n.\n");
	return(FAILURE);
    };
    for(n=0;n<=taylor_n;n++)
	mpfi_c_init(m_sum2[n]);
    mpfi_c_set(m_sum2[0],one);    // m_sum2[0] <- 1

    m_sum1_pos=malloc((taylor_n+1) * sizeof(mpfi_c_t));
    if(m_sum1_pos==NULL)
    {
	printf("Error allocating memory for n_sum1_pos.\n");
	return(FAILURE);
    };
    for(n=0;n<=taylor_n;n++)
	mpfi_c_init(m_sum1_pos[n]);
    mpfi_c_set(m_sum1_pos[0],one);

    m_sum1_neg=malloc((taylor_n+1) * sizeof(mpfi_c_t));
    if(m_sum1_neg==NULL)
    {
	printf("Error allocating memory for n_sum1_neg.\n");
	return(FAILURE);
    };
    for(n=0;n<=taylor_n;n++)
	mpfi_c_init(m_sum1_neg[n]);
    mpfi_c_set(m_sum1_neg[0],one);


   

    mpfi_init(h);
    mpfi_init(tmp);
    mpfi_c_init(s);
    mpfi_c_init(ks);
    mpfi_c_init(zero);
    mpfi_c_init(F_s);
    mpfi_c_init(ks_sqr);
    mpfi_c_init(n_sum_pos);
    mpfi_c_init(n_sum_neg);
    mpfi_c_init(c_tmp1);
//    mpfi_c_init(g_rho);
    mpfi_c_init(g_rho_delta);
    mpfi_c_init(g_rho_delta_sum);
    mpfi_set_ui(tmp,0);
    mpfi_c_set_re(zero,tmp);
    mpfi_c_set_im(zero,tmp);
    mpfi_c_set(g_rho_delta_sum,zero);
    mpfi_c_set_re(s,half);
    mpfi_set_d(s->im,max_height);      // s=1/2+i*max_height

    mpfr_init(h_fr);
    mpfr_init(two_h_max);
    mpfi_init(h_max_i);
    mpfr_add(two_h_max,h_max,h_max,GMP_RNDN);
    mpfi_set_fr(h_max_i,h_max);
    
    mpfi_get_fr(h_fr,s->im);
    while(mpfr_cmp(h_fr,two_h_max)>0)  // s is too far from 1/2 for one leap
    {
        mpfi_sub(s->im,s->im,h_max_i);
	G(g_rho_delta,s,h_max_i);
	mpfi_c_add(g_rho_delta_sum,g_rho_delta_sum,g_rho_delta);
	mpfi_sub(s->im,s->im,h_max_i);
	mpfi_get_fr(h_fr,s->im);
    };

    mpfi_div_ui(s->im,s->im,2);
    G(g_rho_delta,s,s->im);
    mpfi_c_add(g_rho_delta_sum,g_rho_delta_sum,g_rho_delta);
    printf("G(1/2+%f)-G(1/2) is ",max_height);
    mpfi_c_print(g_rho_delta_sum);


// clean everything up and return

    mpfi_clear(half);
    mpfi_clear(h);
    mpfi_clear(tmp);
    mpfi_c_clear(one);
    mpfi_c_clear(s);
    mpfi_c_clear(ks);
    mpfi_c_clear(F_s);
    mpfi_c_clear(ks_sqr);
    mpfi_c_clear(n_sum_pos);
    mpfi_c_clear(n_sum_neg);
    mpfi_c_clear(g_rho_delta);
    mpfi_c_clear(g_rho_delta_sum);
    mpfi_c_clear(c_tmp1);
//    mpfi_c_clear(g_rho);
    mpfr_clear(h_fr);
    mpfr_clear(two_h_max);
    mpfi_clear(h_max_i);
    for(n=0;n<=taylor_n;n++)
    {
	mpz_clear(nbang[n]);
	mpfi_c_clear(ihk_n[n]);
	mpfi_c_clear(m_sum2[n]);
	mpfi_c_clear(m_sum1_neg[n]);
	mpfi_c_clear(m_sum1_pos[n]);
    };
    for(n=taylor_n+1;n<=taylor_n+taylor_n;n++)
    {
	mpz_clear(nbang[n]);
	mpfi_c_clear(ihk_n[n]);
    };
    mpfi_c_clear(ihk_n[taylor_n+taylor_n+1]);
    
    free(nbang);
    free(ihk_n);
    free(ks_2n_kks);
    free(m_sum2);
    free(m_sum1_pos);
    free(m_sum1_neg);

    return(SUCCESS);
};
    
void print_usage()
{
	printf("Usage: G <bits_prec> <x> <lambda> <taylor_n> <h_max>\n%lu <= bits_prec <= %lu\n",
	       MPFR_PREC_MIN,MPFR_PREC_MAX);
	printf("x >= 2\n");
	printf("lambda > 0\n");
	printf("taylor_n > 0\n");
	printf("0<h_max<=%f\n",H_MAX);
};

int main(int argc, char **argv)
{

    int prec;

/*  check all the command line arguments are ok, if not print message
    and exit sharpish */


    if(argc!=6)
    {
	print_usage();
	return(QUIT_CODE);
    };

    prec=atoi(argv[1]);
    if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    {
	print_usage();
	return(QUIT_CODE);
    };
    
/* set the default precision for mpfr, ergo mpfi and mpfi_c as well */ 
    printf("Setting Default Precision to %d bits.\n",prec);
    mpfr_set_default_prec(prec);
    mpfi_c_setup();
    mpfr_init(lambda);
    mpfi_init(lam_sqr);
    mpfi_init(lam_sqr_by_2);
    mpz_init(x);
    mpfi_init(lnx);

    if(mpz_set_str(x,argv[2],10)!=SUCCESS)
    {
	print_usage();
	return(QUIT_CODE);
    };
    if(mpz_cmp_ui(x,2)<0)
    {
	print_usage();
	return(QUIT_CODE);
    };

/*  The value of lambda is rounded to the nearest exact float subject
    the precision specified. Since lambda is a parameter, not
    a result of a calculation, this won't matter so long as we use
    exactly the same lambda in the prime sieve as well.
*/

    if(mpfr_set_str(lambda,argv[3],10,GMP_RNDN)!=SUCCESS)
    {
	print_usage();
	return(QUIT_CODE);
    };
    if (mpfr_sgn(lambda)!=1)
    {
	print_usage();
	return(QUIT_CODE);
    };

    if((taylor_n=atoi(argv[4]))<=0)
    {
	print_usage();
	return(QUIT_CODE);
    };

    mpfr_init(h_max);
    if(mpfr_set_str(h_max,argv[5],10,GMP_RNDN)!=SUCCESS)
      {
	print_usage();
	return(QUIT_CODE);
      };
    if((mpfr_sgn(h_max)<=0)||(mpfr_cmp_d(h_max,H_MAX)>0))
      {
	print_usage();
	return(QUIT_CODE);
      };

    printf("lambda set to ");
    mpfr_print(lambda);
    printf("taylor_n set to %d.\n",taylor_n);
    printf("h_max set to ");
    mpfr_print(h_max);


/* initialise our interval variables */
    printf("Initialising variables.\n");


    

   


    mpfi_set_fr(lam_sqr,lambda);
    mpfi_mul(lam_sqr,lam_sqr,lam_sqr);
    mpfi_div_ui(lam_sqr_by_2,lam_sqr,2);
    printf("Lambda^2/2 set to ");
    mpfi_print(lam_sqr_by_2);
    mpfi_set_z(lnx,x);
    mpfi_log(lnx,lnx);
    printf("ln(x) set to  ");
    mpfi_print(lnx);

	if(g_half()!=SUCCESS)
	   return(QUIT_CODE);

/* not bothering to clear data structures on exit */    
    return(SUCCESS);
}

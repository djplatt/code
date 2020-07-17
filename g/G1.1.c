/*

File: G1.1.c

Created: 4 March 2008

Version: 1.1

Last Modified: 5 March 2008

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

Build instructions: gcc -oG1.1 G1.1.c -O2 -lmpfi -lmpfr -lgmp

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */


/*

Calculates an estimate of the sum of -G(rho) where rho is a non-trivial zero
of the Riemann Zeta Function, plus G(1) etc. etc. */

#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "mpfi_c.c"

#define SUCCESS 0
#define QUIT_CODE 0
#define FAILURE 1
#define LOAD_ZEROS_ERROR 1
#define ZEROS_FILE "zeros1.dat"
/* the first line of this file is an arbitary precision scale factor (a power
   of 10) by which each of the non-trivial zeros has been multiplied to make
   it an integer. The next line is an unsigned long integer tollerance of the
   least significant digit (plus and minus). The next line is an unsigned
   long long integer, how many zeros the file contains. Then each line is the
   imaginary part of a non-trivial zero above the real axis, in order of
   increasing Im(rho). */



#define TRUE (1==1)
#define FALSE (0==1)
#define H_MAX 0.5       /* absolute upper limit for h_max otherwise taylor won't converge */




/* global variables for commonly used values initialised in main */

mpfr_t lambda;
mpfi_t lam_sqr,lam_sqr_by_2,lnx;
mpz_t x;
int taylor_n;
mpfi_t h_max,two_h_max;
mpfi_t g_one_height;


/* the imaginary parts of the zeros will be kept in here */

mpfi_t *zeros;

/* and this is how many we have */
long long unsigned int NUM_ZEROS;
/* both initialised by load_zeros */

int load_zeros()
{
    int i;
    mpz_t rho_up,rho_down,scale;
    unsigned long int tol;              // the tollerance
    FILE *zero_file;

    printf("About to Read %s.\n",ZEROS_FILE);

    zero_file=fopen(ZEROS_FILE,"r");
    if (zero_file==NULL)
    {
	printf("Error opening %d\n",ZEROS_FILE);
	return(FAILURE);
    };
    
    mpz_init(scale);
    mpz_inp_str(scale,zero_file,10);

    fscanf(zero_file,"%lu%Lu",&tol,&NUM_ZEROS);

    zeros=malloc(NUM_ZEROS * sizeof(mpfi_t));
    if (zeros==NULL)
    {
        printf("Error allocating memory for zeros vector.\n");
	return(FAILURE);
    };

    mpz_init(rho_up);
    mpz_init(rho_down);

    for(i=0;i<NUM_ZEROS;i++)
    {
	mpz_inp_str(rho_up,zero_file,10);
	mpz_add_ui(rho_up,rho_up,tol);
	mpz_sub_ui(rho_down,rho_up,tol+tol);
	mpfi_init(zeros[i]);
	mpfi_interv_z(zeros[i],rho_down,rho_up);
	mpfi_div_z(zeros[i],zeros[i],scale);
    };

    fclose(zero_file);

//    mpz_clear(rho_up);
//    mpz_clear(rho_down);
//    mpz_clear(scale);
    printf("%Ld zeros read.\n",NUM_ZEROS);
    printf("Last zero at height ");
    mpfi_print(zeros[NUM_ZEROS-1]);

    return(SUCCESS);
};

void free_zeros()
{
    long long unsigned int i;
/*
    for(i=0;i<NUM_ZEROS;i++)
        mpfi_clear(zeros[i]);

    free(zeros);
*/
};


void F (mpfi_c_ptr res,mpfi_c_ptr s)
{
    mpfi_c_mul(F_tmp1,s,s);
    mpfi_c_mul_i(F_tmp1,F_tmp1,lam_sqr_by_2);
    mpfi_c_mul_i(F_tmp2,s,lnx);
    mpfi_c_add(F_tmp1,F_tmp1,F_tmp2);
    mpfi_c_exp(F_tmp1,F_tmp1);
    mpfi_c_div(res,F_tmp1,s);
};


/* global variable used by G
   initialised in init_vectors()
   destroyed in free_vectors()
*/

mpfi_c_t ks,F_s,n_sum_pos;
mpfi_c_t n_sum_neg,ks_sqr,c_tmp1,zero;
mpz_t *nbang;
mpfi_c_t *ihk_n;
mpfi_c_t *ks_2n_kks;        // k^(2n+2)s^(2n+1)
mpfi_c_t *m_sum2;    // sigma (lambda^2*s^2/2)^n/n!
mpfi_c_t *m_sum1_pos;
mpfi_c_t *m_sum1_neg;

/* this sets err to be the lower and upper bound for the real integral error,
   i.e. the error introduced by taking G(1+g_one_height*i)=0
   
*/

void calc_one_integral_error(mpfi_ptr err)
{
    mpfi_t ans,tmp1;

    mpfi_init(ans);
    mpfi_init(tmp1);

    mpfi_sqr(tmp1,g_one_height);
    mpfi_mul(ans,tmp1,lam_sqr_by_2);
    mpfi_neg(ans,ans);
    mpfi_exp(ans,ans);
//    printf("exp(-(T*lam)^2/2 is ");
//    mpfi_print(ans);
    mpfi_div(ans,ans,lam_sqr_by_2);
    mpfi_div(ans,ans,tmp1);

    mpfi_exp(tmp1,lam_sqr_by_2);
    mpfi_mul_z(tmp1,tmp1,x);
    mpfi_mul(ans,ans,tmp1);
    mpfr_neg(&ans->left,&ans->left,GMP_RNDD);    

    mpfi_add(err,err,ans);

    printf("Error truncating integral on Re(s)=1 line was ");
    mpfi_print(ans);

//    mpfi_clear(ans);
//    mpfi_clear(tmp1);
    
};




/* this sets err to be the lower and upper bound for the real integral error,
   i.e. the error introduced by taking G(1/2+rho_max*i)=0
   
*/
void calc_half_integral_error(mpfi_ptr err)
{
    mpfi_t ans,tmp1;


    mpfi_init(ans);
    mpfi_init(tmp1);

    mpfi_sqr(tmp1,zeros[NUM_ZEROS-1]);
    mpfi_mul(ans,tmp1,lam_sqr_by_2);
    mpfi_neg(ans,ans);
    mpfi_exp(ans,ans);
    mpfi_div(ans,ans,lam_sqr_by_2);
    mpfi_div(ans,ans,tmp1);

    mpfi_div_ui(tmp1,lam_sqr_by_2,4);
    mpfi_exp(tmp1,tmp1);
    mpfi_mul(ans,ans,tmp1);
    
    mpfi_set_z(tmp1,x);
    mpfi_sqrt(tmp1,tmp1);
    mpfi_mul(ans,ans,tmp1);
    mpfr_neg(&ans->left,&ans->left,GMP_RNDD);    

    mpfi_add(err,err,ans);

    printf("Error truncating integral on Re(s)=1/2 line was ");
    mpfi_print(ans);

//    mpfi_clear(ans);
//    mpfi_clear(tmp1);
    
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

/* adds the error due to the remaining integral at Re(s)=-1 */
void re_minus_1_error(mpfi_ptr err)
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
	};

	mpfi_c_exp(ans,ihk_n[1]);
	mpfi_c_mul(ans,ans,n_sum_pos);
	mpfi_c_neg(c_tmp1,ihk_n[1]);
	mpfi_c_exp(c_tmp1,c_tmp1);
	mpfi_c_mul(c_tmp1,c_tmp1,n_sum_neg);
	mpfi_c_sub(ans,ans,c_tmp1);
	F(F_s,s);
	mpfi_c_mul(ans,F_s,ans);
	calc_taylor_error(ans,s,h);


}; /* G */

int init_vectors()
{
    mpfi_c_t one;
    int n;

    mpfi_c_init(one);

    mpfi_c_init(ks);
    mpfi_c_init(F_s);
    mpfi_c_init(ks_sqr);
    mpfi_c_init(n_sum_pos);
    mpfi_c_init(n_sum_neg);
    mpfi_c_init(c_tmp1);
    mpfi_c_init(zero);


    
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
    ihk_n=malloc((taylor_n+taylor_n+2) * sizeof(mpfi_c_t));
    if(ihk_n==NULL)
    {
	printf("Error allocating memory for ihk_n.\n");
	return(FAILURE);
    };
    for(n=0;n<taylor_n+taylor_n+2;n++)
	mpfi_c_init(ihk_n[n]);

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


//    mpfi_c_clear(one);

    return(SUCCESS);
};

void free_vectors()
{
    int n;
/*
    mpfi_c_clear(zero);
    mpfi_c_clear(ks);
    mpfi_c_clear(F_s);
    mpfi_c_clear(n_sum_pos);
    mpfi_c_clear(n_sum_neg);
    mpfi_c_clear(ks_sqr);
    mpfi_c_clear(c_tmp1);

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
*/
};


int g_sum(mpfi_ptr g_rho_sum)
{
    unsigned long long int i;
    unsigned long int n,m;
    mpfi_t half,h,upper,lower;
    mpfi_c_t s,g_rho,g_rho_delta;
    mpfi_c_t g_rho_delta_sum;

    mpfi_init(half);
    mpfi_init(h);
    mpfi_init(upper);
    mpfi_init(lower);
    mpfi_c_init(s);
    mpfi_c_init(g_rho);
    mpfi_c_init(g_rho_delta);
    mpfi_c_init(g_rho_delta_sum);

/* assign 1/2 to half */
    if(mpfi_set_d(half,0.5)!=MPFI_FLAGS_BOTH_ENDPOINTS_EXACT)
    {
	printf("Inconceivable error setting up 1/2.\n");
	return(FAILURE);
    };


    if(load_zeros()!=SUCCESS)
	return(FAILURE);

    mpfi_c_set(g_rho,zero);
    mpfi_set_ui(g_rho_sum,0);
    mpfi_c_set_re(s,half);

    for(i=NUM_ZEROS-1;i>0;i--)
/*  looking at the integral between zeros[i] and zeros[i-1]
*/
    {
	mpfi_c_set(g_rho_delta_sum,zero);
	mpfi_set(upper,zeros[i]);
	mpfi_sub(s->im,upper,h_max);
	mpfi_sub(lower,s->im,h_max);
	while(mpfi_cmp(lower,zeros[i-1])>0)
	{
	    G(g_rho_delta,s,h_max);
	    mpfi_c_add(g_rho_delta_sum,g_rho_delta_sum,g_rho_delta);
	    mpfi_sub(upper,upper,two_h_max);
	    mpfi_sub(s->im,upper,h_max);
	    mpfi_sub(lower,s->im,h_max);
	};
	mpfi_sub(h,upper,zeros[i-1]);
	mpfi_div_ui(h,h,2);
	mpfi_sub(s->im,upper,h);
	G(g_rho_delta,s,h);
	mpfi_c_add(g_rho,g_rho_delta_sum,g_rho_delta);
	mpfi_sub(g_rho_sum,g_rho_sum,g_rho->re);
	mpfi_sub(g_rho_sum,g_rho_sum,g_rho->re);
    };

    free_zeros();
/*
    mpfi_clear(half);
    mpfi_clear(h);
    mpfi_clear(upper);
    mpfi_clear(lower);
    mpfi_c_clear(s);
    mpfi_c_clear(g_rho_delta);
    mpfi_c_clear(g_rho_delta_sum);
    mpfi_c_clear(g_rho);
*/
    return(SUCCESS);
};

int g_one(mpfi_ptr g1)
{

    mpfi_c_t s,g_one_delta,g_one;
    mpfi_t upper,lower,h;


    mpfi_c_init(s);
    mpfi_c_init(g_one_delta);
    mpfi_c_init(g_one);


    mpfi_init(upper);
    mpfi_init(lower);
    mpfi_init(h);

    mpfi_set_ui(s->re,1);
    mpfi_c_set(g_one,zero);
    mpfi_set(upper,g_one_height);
    mpfi_sub(s->im,upper,h_max);
    mpfi_sub(lower,s->im,h_max);
    while(mpfi_cmp_ui(lower,0)>0)
    {
	G(g_one_delta,s,h_max);
	mpfi_c_add(g_one,g_one,g_one_delta);
	mpfi_sub(upper,upper,two_h_max);
	mpfi_sub(s->im,upper,h_max);
	mpfi_sub(lower,s->im,h_max);
    };
    mpfi_div_ui(h,upper,2);
    mpfi_set(s->im,h);
    G(g_one_delta,s,h);
    mpfi_c_add(g_one,g_one,g_one_delta);
    mpfi_neg(g1,g_one->re);

//    mpfi_c_clear(s);
//    mpfi_c_clear(g_one_delta);
//    mpfi_c_clear(g_one);

//    mpfi_clear(upper);
//   mpfi_clear(lower);
//    mpfi_clear(h);


    return(SUCCESS);
}
    
void print_usage()
{
	printf("Usage: G <bits_prec> <x> <lambda> <taylor_n> <h_max> <G(1) height>\n%lu <= bits_prec <= %lu\n",
	       MPFR_PREC_MIN,MPFR_PREC_MAX);
	printf("x >= 2\n");
	printf("lambda > 0\n");
	printf("taylor_n > 0\n");
	printf("0<h_max<=%f\n",H_MAX);
	printf("G(1) height >0\n");
};

int main(int argc, char **argv)
{

    int prec;
    mpfi_t g1,ln2,err_term;
    mpfi_t g_rho_sum;
    mpfr_t h_max_fr;

/*  g_rho_sum will contain the sum of all the G(rho) terms, plus the
    error because we have ignored zeros > NUM_ZEROS, plus the error
    because we assume G(biggest rho) is zero.

    g1 will contain our estimate of G(1)

    ln2 will contain natural log(2)

    err_term is the bounds for the integral of F(s) with Im(s)=-1

    if we add all the lower bounds and all the upper bounds of these real
    intervals, we get the interval containing our result. */

/*  check all the command line arguments are ok, if not print message
    and exit sharpish */


    if(argc!=7)
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
    
    mpfi_c_setup(prec);

    mpz_init(x);
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

    mpfr_init(lambda);

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

    mpfr_init(h_max_fr);
    if(mpfr_set_str(h_max_fr,argv[5],10,GMP_RNDN)!=SUCCESS)
      {
	print_usage();
	return(QUIT_CODE);
      };
    if((mpfr_sgn(h_max_fr)<=0)||(mpfr_cmp_d(h_max_fr,H_MAX)>0))
      {
	print_usage();
	return(QUIT_CODE);
      };

    mpfi_init(h_max);
    mpfi_init(two_h_max);
    mpfi_set_fr(h_max,h_max_fr);
    mpfi_add(two_h_max,h_max,h_max);

    if(mpfr_set_str(h_max_fr,argv[6],10,GMP_RNDN)!=SUCCESS)
	{
	    print_usage();
	    return(QUIT_CODE);
	};
    if(mpfr_sgn(h_max_fr)<=0)
	{
	    print_usage();
	    return(QUIT_CODE);
	};
    mpfi_init(g_one_height);
    mpfi_set_fr(g_one_height,h_max_fr);

    printf("lambda set to ");
    mpfr_print(lambda);
    printf("taylor_n set to %d.\n",taylor_n);
    printf("h_max set to ");
    mpfi_print(h_max);
    printf("G(1) height set to ");
    mpfi_print(g_one_height);


/* initialise our interval variables */
    printf("Initialising variables.\n");
    mpfi_init(g_rho_sum);
    mpfi_init(g1);
    mpfi_init(ln2);
    mpfi_init(err_term);
    mpfi_init(lam_sqr);
    mpfi_init(lam_sqr_by_2);
    mpfi_init(lnx);


    mpfi_set_fr(lam_sqr,lambda);
    mpfi_sqr(lam_sqr,lam_sqr);
    mpfi_div_ui(lam_sqr_by_2,lam_sqr,2);
    printf("Lambda^2/2 set to ");
    mpfi_print(lam_sqr_by_2);
    mpfi_set_z(lnx,x);
    mpfi_log(lnx,lnx);
    printf("ln(x) set to  ");
    mpfi_print(lnx);


    
/* call standard mpfi routine to assign ln2:=ln(2) */
    mpfi_const_log2(ln2);
    printf("ln(2) set to ");
    mpfi_print(ln2);

    printf("Initialising vectors.\n");
    if(init_vectors()!=SUCCESS)
	return(QUIT_CODE);

    printf("Entering g_sum.\n");
    if(g_sum(g_rho_sum)!=SUCCESS)
	return(QUIT_CODE);

    printf("Sum of G(rho) is ");mpfi_print(g_rho_sum);

    printf("Entering g_one.\n");
    if(g_one(g1)!=SUCCESS)
	return(QUIT_CODE);



    printf("G(1) is ");mpfi_print(g1);

    

    mpfi_sub(g1,g1,g_rho_sum);
    mpfi_sub(g1,g1,ln2);
    re_minus_1_error(g1);
    calc_one_integral_error(g1);
    calc_half_integral_error(g1);
    calc_trunc_error(g1);

    
    printf("I1 is in the interval ");
    mpfi_print(g1);

    free_vectors();

/* not bothering to clear data structures on exit */    
    return(SUCCESS);
}

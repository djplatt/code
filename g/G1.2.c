/*

File: G1.2.c

Created: 4 March 2008

Version: 1.2

Last Modified: 2 April 2008

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

V 1.0 Initial implementation
V 1.1 Control over maximum size of h added
      Calculation of G(1) added
      mpfi_c routines moved to mpfi_c.c
V 1.2 Revised Taylor Expansion used and taylor error added

Build instructions: gcc -oG1.2 G1.2.c -O2 -lmpfi -lmpfr -lgmp

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
mpfi_t g_one_height,taylor_err;


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
    printf("Last zero at height\n");
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

mpfi_c_t F_s,zero,*ihk_n,k,k_n,one,*ihk_n_sum_pos,*ihk_n_sum_neg;
mpfi_c_t *s_minus_n,mit,msum,nsum_pos,nsum_neg,nsum;
mpfi_t *lam_sqr_by_2_n,t_tmp1,t_tmp2,t_tmp5;
mpfi_c_t t_tmp3,t_tmp4;
mpz_t *nbang;



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

    printf("Error truncating integral on Re(s)=1 line was\n");
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

// since we estimate g(rho) 2*NUM_ZEROS times,
// that error accumulates

    mpfi_mul_ui(ans,ans,NUM_ZEROS+NUM_ZEROS);


    mpfi_add(err,err,ans);

    printf("Error truncating integral on Re(s)=1/2 line was\n");
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
   and s,-h to taylor_err
   
*/

/* adds the real taylor error at s and h,-h to err */
void calc_taylor_error(mpfi_c_ptr s,mpfi_ptr h)
{
    mpfi_c_exp(t_tmp3,ihk_n[1]);
    mpfi_c_neg(t_tmp4,ihk_n[1]);
    mpfi_c_exp(t_tmp4,t_tmp4);
    mpfi_c_sub(t_tmp3,t_tmp3,t_tmp4);
    mpfi_c_mul(t_tmp3,t_tmp3,F_s);
    mpfi_c_abs(t_tmp5,t_tmp3);
    mpfi_c_abs(t_tmp1,s);
    mpfi_sub(t_tmp2,t_tmp1,h);
    mpfi_div(t_tmp1,h,t_tmp1);
    mpfi_log(t_tmp1,t_tmp1);
    mpfi_mul_ui(t_tmp1,t_tmp1,taylor_n+1);
    mpfi_exp(t_tmp1,t_tmp1);
    mpfi_mul(t_tmp1,t_tmp1,h);
    mpfi_div(t_tmp1,t_tmp1,t_tmp2);
    mpfi_mul(t_tmp1,t_tmp1,t_tmp5);

    mpfr_neg(&t_tmp1->left,&t_tmp1->left,GMP_RNDD);
    mpfi_add(taylor_err,taylor_err,t_tmp1);
	
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
  	unsigned int n,m;

	F(F_s,s);	// calculate F(s)


// set ihk_n[1]=i*h*(lambda^2*s+ln(x))
// set ihk_n_sum_neg[n]=sum 0..n ihk_n[n]
// set ihk_n_sum_pos[n]=sum 0..n (-1)^n*ihk_n[n]
	mpfi_c_set(k,s);
	mpfi_c_mul_i(k,k,lam_sqr);
	mpfi_c_add_re(k,k,lnx);
	mpfi_c_mul_i(ihk_n[1],k,h);
	mpfi_c_muli(ihk_n[1]);

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

	mpfi_c_neg(msum,ihk_n[1]);       // msum = -ihk
	mpfi_c_exp(msum,msum);		 // msum = exp(-ihk)
	mpfi_c_mul(msum,msum,nsum_neg);  // msum = exp(-ihk)*sum
	mpfi_c_exp(ans,ihk_n[1]);        // ans  = exp(ihk)
	
	mpfi_c_mul(ans,ans,nsum_pos);    // ans  = exp(ihk)*sum
	
	mpfi_c_sub(ans,msum,ans);
	
	mpfi_c_mul(ans,ans,F_s);



	calc_taylor_error(s,h);


}; /* G */

int init_vectors()
{
    
    int n;

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
    mpfi_init(t_tmp1);
    mpfi_init(t_tmp2);
    mpfi_init(t_tmp5);
    mpfi_c_init(t_tmp3);
    mpfi_c_init(t_tmp4);
    
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

void free_vectors()
{
    int n;

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
	mpfi_c_set(g_rho_delta_sum,g_rho);
	mpfi_set(upper,zeros[i]);
	mpfi_sub(s->im,upper,h_max);
	mpfi_sub(lower,s->im,h_max);
	while(mpfi_cmp(lower,zeros[i-1])>0)
	{
	    G(g_rho_delta,s,h_max);
//	    printf("s was ");mpfi_c_print(s);
//          printf("g_rho_delta is ");mpfi_c_print(g_rho_delta);
	    mpfi_c_add(g_rho_delta_sum,g_rho_delta_sum,g_rho_delta);
//	    printf("g_rho_delta_sum is ");mpfi_c_print(g_rho_delta_sum);
	    mpfi_sub(upper,upper,two_h_max);
	    mpfi_sub(s->im,upper,h_max);
	    mpfi_sub(lower,s->im,h_max);
	};
	mpfi_sub(h,upper,zeros[i-1]);
	mpfi_div_ui(h,h,2);
	mpfi_sub(s->im,upper,h);
	G(g_rho_delta,s,h);
//	printf("s was ");mpfi_c_print(s);
//    printf("g_rho_delta is ");mpfi_c_print(g_rho_delta);
	mpfi_c_add(g_rho,g_rho_delta_sum,g_rho_delta);
//	printf("G(");mpfi_out_str(stdout,10,0,zeros[i-1]);printf(") was\n");mpfi_c_print(g_rho);
	mpfi_add(g_rho_sum,g_rho_sum,g_rho->re);
	mpfi_add(g_rho_sum,g_rho_sum,g_rho->re);
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

/* Calculate an estimate of Re(G(1)) */

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
    mpfi_set(g1,g_one->re);
    printf("Imaginary part of G(1) was ");
    mpfi_print(g_one->im);
    printf("Real part of G(1) was ");
    mpfi_print(g_one->re);
//    mpfi_c_clear(s);
//    mpfi_c_clear(g_one_delta);
//    mpfi_c_clear(g_one);

//    mpfi_clear(upper);
//   mpfi_clear(lower);
//    mpfi_clear(h);


    return(SUCCESS);
};
    
void print_usage()
{
	printf("Usage: G <bits_prec> <x> <lambda> <taylor_n> <h_max> <G(1) height>\n%lu <= bits_prec <= %lu\n",
	       MPFR_PREC_MIN,MPFR_PREC_MAX);
	printf("x >= 2\n");
	printf("lambda > 0\n");
	printf("taylor_n > 0, taylor_n even\n");
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

    if(((taylor_n=atoi(argv[4]))<=0)||(taylor_n%2==1))
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

    printf("lambda set to\n");
    mpfr_print(lambda);
    printf("taylor_n set to %d.\n",taylor_n);
    printf("h_max set to\n");
    mpfi_print(h_max);
    printf("G(1) height set to\n");
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
    mpfi_init(taylor_err);


    mpfi_set_fr(lam_sqr,lambda);
    mpfi_sqr(lam_sqr,lam_sqr);
    mpfi_div_ui(lam_sqr_by_2,lam_sqr,2);
    printf("Lambda^2/2 set to\n");
    mpfi_print(lam_sqr_by_2);
    mpfi_set_z(lnx,x);
    mpfi_log(lnx,lnx);
    printf("ln(x) set to\n");
    mpfi_print(lnx);
    mpfi_set_ui(taylor_err,0);


    
/* call standard mpfi routine to assign ln2:=ln(2) */
    mpfi_const_log2(ln2);
    printf("ln(2) set to\n");
    mpfi_print(ln2);

    printf("Initialising vectors.\n");
    if(init_vectors()!=SUCCESS)
	return(QUIT_CODE);

    printf("Entering g_sum.\n");
    if(g_sum(g_rho_sum)!=SUCCESS)
	return(QUIT_CODE);


    printf("Sum of G(rho) is\n ");
    mpfi_print(g_rho_sum);
    printf("Taylor error calculating sum G(rho)\n ");
    mpfi_print(taylor_err);
    mpfi_add(g_rho_sum,g_rho_sum,taylor_err);
    mpfi_set_ui(taylor_err,0);

    printf("Entering g_one.\n");
    if(g_one(g1)!=SUCCESS)
	return(QUIT_CODE);



    printf("G(1) is\n");
    mpfi_print(g1);
    printf("Taylor error calculating G(1)\n ");
    mpfi_print(taylor_err);
    mpfi_add(g1,g1,taylor_err);    

    mpfi_sub(g1,g1,g_rho_sum);
    mpfi_sub(g1,g1,ln2);
    re_minus_1_error(g1);
    calc_one_integral_error(g1);
    calc_half_integral_error(g1);
    calc_trunc_error(g1);

    
    printf("I1 is in the interval\n");
    mpfi_print(g1);

    free_vectors();

/* not bothering to clear data structures on exit */    
    return(SUCCESS);
}

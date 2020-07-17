/*

File: G.c

Created: 25 February 2008

Version: 1.0

Last Modified: 29 February 2008

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

Build instructions: gcc -oG G.c -O2 -lmpfi -lmpfr -lgmp

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

#define SUCCESS 0
#define QUIT_CODE 0
#define FAILURE 1
#define LOAD_ZEROS_ERROR 1
#define ZEROS_FILE "zeros1.dat"
#define TRUE (1==1)
#define FALSE (0==1)

/* the first line of this file is an arbitary precision scale factor (a power
   of 10) by which each of the non-trivial zeros has been multiplied to make
   it an integer. The next line is an unsigned long integer tollerance of the
   least significant digit (plus and minus). The next line is an unsigned
   long long integer, how many zeros the file contains. Then each line is the
   imaginary part of a non-trivial zero above the real axis, in order of
   increasing Im(rho). */


/* define a complex version of an MPFI interval */
typedef struct{
    mpfi_t re;
    mpfi_t im;
} _mpfi_c_struct;

typedef _mpfi_c_struct mpfi_c_t[1];
typedef _mpfi_c_struct *mpfi_c_ptr;

/* initialisation */
void mpfi_c_init(mpfi_c_ptr z)
{
    mpfi_init(z->re);
    mpfi_init(z->im);
};

/* clearing */
void mpfi_c_clear(mpfi_c_ptr z)
{
    mpfi_clear(z->re);
    mpfi_clear(z->im);
};

/* swapping */

void mpfi_c_swap(mpfi_c_ptr z1,mpfi_c_ptr z2)
{
    mpfi_swap(z1->re,z2->re);
    mpfi_swap(z1->im,z2->im);
};

/* print one */
void mpfi_c_print(mpfi_c_ptr z)
{
    mpfi_out_str(stdout,10,0,z->re);
    printf(" + ");
    mpfi_out_str(stdout,10,0,z->im);
    printf("i\n");
};

mpfi_print(mpfi_ptr x)
{
    mpfi_out_str(stdout,10,0,x);
    printf("\n");
};

mpfr_print(mpfr_ptr x)
{
    mpfr_out_str(stdout,10,0,x,GMP_RNDN);
    printf("\n");
};

void mpfi_c_set(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_set(z1->re,z2->re);
    mpfi_set(z1->im,z2->im);
};

void mpfi_c_set_re(mpfi_c_ptr z,mpfi_ptr x)
{
    mpfi_set(z->re,x);
};

void mpfi_c_set_re_d(mpfi_c_ptr z,double x)
{
    mpfi_set_d(z->re,x);
};

void mpfi_c_set_im(mpfi_c_ptr z, mpfi_ptr x)
{
    mpfi_set(z->im,x);
};



void mpfi_c_add(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_add(z1->re,z2->re,z3->re);
    mpfi_add(z1->im,z2->im,z3->im);
};

void mpfi_c_sub(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_sub(z1->re,z2->re,z3->re);
    mpfi_sub(z1->im,z2->im,z3->im);
};

void mpfi_c_add_re(mpfi_c_ptr z1,mpfi_c_ptr z2,mpfi_ptr x)
{
    mpfi_add(z1->re,z2->re,x);
    mpfi_set(z1->im,z2->im);
};

void mpfi_c_add_ui(mpfi_c_ptr z1,mpfi_c_ptr z2, unsigned long int n)
{
    mpfi_add_ui(z1->re,z2->re,n);
    mpfi_set(z1->im,z2->im);
};

/* global varibale for commonly used values */

mpfr_t lambda;
mpfi_t lam_sqr,lam_sqr_by_2,lnx;
mpz_t x;
int taylor_n;


/* temporary variables for mul/div/exp/F
   declared globally to avoid having to initialise
   each time mul etc. gets called */

mpfi_c_t d_tmp3,F_tmp1,F_tmp2;
mpfi_t m_tmp1,m_tmp2,m_tmp3,d_tmp1,d_tmp2,e_tmp;   

void mpfi_c_setup(int prec)
{
/* set the default precision for mprf, ergo mpfi and mpfi_c as well */ 
    printf("Setting Default Precision to %d bits.\n",prec);
    mpfr_set_default_prec(prec);
/* now initialise all global variables.
   we don't bother clearing them down afterwards, just exit. */
    mpfi_init(m_tmp1);
    mpfi_init(m_tmp2);
    mpfi_init(m_tmp3);
    mpfi_init(d_tmp1);
    mpfi_init(d_tmp2);
    mpfi_c_init(d_tmp3);
    mpfi_init(e_tmp);
    mpfi_c_init(F_tmp1);
    mpfi_c_init(F_tmp2);
    mpfr_init(lambda);
    mpfi_init(lam_sqr);
    mpfi_init(lam_sqr_by_2);
    mpz_init(x);
    mpfi_init(lnx);

};

void mpfi_c_mul(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
/*    printf("mpfi_c_mul called with ");
    mpfi_c_print(z1);
    mpfi_c_print(z2);
    mpfi_c_print(z3); */
    mpfi_mul(m_tmp1,z2->re,z3->re);
    mpfi_mul(m_tmp2,z2->im,z3->im);
    mpfi_sub(m_tmp3,m_tmp1,m_tmp2);
    mpfi_mul(m_tmp1,z2->re,z3->im);
    mpfi_mul(m_tmp2,z2->im,z3->re);
    mpfi_add(z1->im,m_tmp1,m_tmp2);
    mpfi_swap(z1->re,m_tmp3);
};

/* multiply (mpfi_c_t z2) by (mpfi_t x) into (mpfi_c_t z1) safely */
void mpfi_c_mul_i(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_ptr x)
{
    mpfi_mul(z1->re,z2->re,x);
    mpfi_mul(z1->im,z2->im,x);
};

void mpfi_c_muli(mpfi_c_ptr z)   /* multiply by i in situ */
{
    mpfi_swap(z->re,z->im);
    mpfi_neg(z->re,z->re);
};

/* multiply (mpfi_c_t z2) by (mpfr_t x) into (mpfi_c_t z1) safely */
void mpfi_c_mul_fr(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfr_ptr x)
{
    mpfi_mul_fr(z1->re,z2->re,x);
    mpfi_mul_fr(z1->im,z2->im,x);
};

void mpfi_c_mul_z (mpfi_c_ptr z1,mpfi_c_ptr z2, mpz_ptr x)
{
    mpfi_mul_z(z1->re,z2->re,x);
    mpfi_mul_z(z1->im,z2->im,x);
};
    

/* z1<-conjugate(z2) safely */
void mpfi_c_conj(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_c_set(z1,z2);
    mpfi_neg(z1->im,z1->im);
};

void mpfi_c_neg(mpfi_c_ptr z1,mpfi_c_ptr z2)
{
    mpfi_neg(z1->im,z2->im);
    mpfi_neg(z1->re,z2->re);
};

void mpfi_c_div(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_sqr(d_tmp1,z3->re);
    mpfi_sqr(d_tmp2,z3->im);
    mpfi_add(d_tmp1,d_tmp1,d_tmp2);
    mpfi_c_conj(d_tmp3,z3);
    mpfi_c_mul(z1,d_tmp3,z2);
    mpfi_div(z1->re,z1->re,d_tmp1);
    mpfi_div(z1->im,z1->im,d_tmp1);
};

void mpfi_c_div_z (mpfi_c_ptr z1,mpfi_c_ptr z2, mpz_ptr x)
{
    mpfi_div_z(z1->re,z2->re,x);
    mpfi_div_z(z1->im,z2->im,x);
};

void mpfi_c_exp(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_exp(e_tmp,z2->re);
    mpfi_cos(z1->re,z2->im);
    mpfi_mul(z1->re,z1->re,e_tmp);
    mpfi_sin(z1->im,z2->im);
    mpfi_mul(z1->im,z1->im,e_tmp);
};


/* the imaginary parts of the zeros will be kept in here */

mpfi_t *zeros;

/* and this is how many we have */
long long unsigned int NUM_ZEROS;


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

    mpz_clear(rho_up);
    mpz_clear(rho_down);
    mpz_clear(scale);

    return(SUCCESS);
};


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


/* this sets err to be the lower and upper bound for the real integral error,
   i.e. the error introduced by taking G(rho_max)=0
   not yet implemented so adds 0+0i (nop)
*/
void calc_integral_error(mpfi_c_ptr err)
{
    
};

/* this adds to err the lower and upper bounds for the real error introduced
   by stopping the sum over rho of G at rho_max
   not yet implemented so adds [0,0] (nop)
*/
void calc_trunc_error(mpfi_c_ptr err)
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

int g_sum(mpfi_c_ptr g_rho_sum)
{
    unsigned long long int i;
    unsigned long int n,m;
    mpfi_t half,h,tmp;
    mpfi_c_t s,ks,zero,F_s,n_sum_pos;
    mpfi_c_t n_sum_neg,ks_sqr,c_tmp,c_tmp1;
    mpfi_c_t g_rho;
    mpz_t *nbang;
    mpfi_c_t *ihk_n;
    mpfi_c_t *ks_2n_kks;        // k^(2n+2)s^(2n+1)
    mpfi_c_t *m_sum2;    // sigma (lambda^2*s^2/2)^n/n!
    mpfi_c_t *m_sum1_pos;
    mpfi_c_t *m_sum1_neg;

    mpfi_init(half);


/* assign 1/2 to half */
    if(mpfi_set_d(half,0.5)!=MPFI_FLAGS_BOTH_ENDPOINTS_EXACT)
    {
	printf("Inconceivable error setting up 1/2.\n");
	return(FAILURE);
    };


    if(load_zeros()!=SUCCESS)
	return(FAILURE);

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
/*    for(n=0;n<=taylor_n+taylor_n;n++)
	{
	    printf("%d! is ",n);mpz_out_str(stdout,10,nbang[n]);printf("\n");
	};
*/
    ihk_n=malloc((taylor_n+taylor_n+2) * sizeof(mpfi_c_t));
    if(ihk_n==NULL)
    {
	printf("Error allocating memory for ihk_n.\n");
	return(FAILURE);
    };
    for(n=0;n<taylor_n+taylor_n+2;n++)
	mpfi_c_init(ihk_n[n]);
    mpfi_set_ui(ihk_n[0]->re,1);
    mpfi_set_ui(ihk_n[0]->im,0);    // ihk_n[0] <- 1+0i=(ihk)^0

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
    mpfi_set_ui(m_sum2[0]->re,1);
    mpfi_set_ui(m_sum2[0]->im,0);    // m_sum2[0] <- 1

    m_sum1_pos=malloc((taylor_n+1) * sizeof(mpfi_c_t));
    if(m_sum1_pos==NULL)
    {
	printf("Error allocating memory for n_sum1_pos.\n");
	return(FAILURE);
    };
    for(n=0;n<=taylor_n;n++)
	mpfi_c_init(m_sum1_pos[n]);
    mpfi_set_ui(m_sum1_pos[0]->re,1);
    mpfi_set_ui(m_sum1_pos[0]->im,0);

    m_sum1_neg=malloc((taylor_n+1) * sizeof(mpfi_c_t));
    if(m_sum1_neg==NULL)
    {
	printf("Error allocating memory for n_sum1_neg.\n");
	return(FAILURE);
    };
    for(n=0;n<=taylor_n;n++)
	mpfi_c_init(m_sum1_neg[n]);
    mpfi_set_ui(m_sum1_neg[0]->re,1);
    mpfi_set_ui(m_sum1_neg[0]->im,0);

    calc_integral_error(g_rho_sum);
    calc_trunc_error(g_rho_sum);

    mpfi_init(h);
    mpfi_init(tmp);
    mpfi_c_init(s);
    mpfi_c_init(ks);
    mpfi_c_init(zero);
    mpfi_c_init(F_s);
    mpfi_c_init(ks_sqr);
    mpfi_c_init(n_sum_pos);
    mpfi_c_init(n_sum_neg);
    mpfi_c_init(c_tmp);
    mpfi_c_init(c_tmp1);
    mpfi_c_init(g_rho);
    mpfi_set_ui(tmp,0);
    mpfi_c_set_re(zero,tmp);
    mpfi_c_set_im(zero,tmp);
    mpfi_c_set(g_rho,zero);
    mpfi_c_set(g_rho_sum,zero);
    mpfi_c_set_re(s,half);

    for(i=NUM_ZEROS-1;i>0;i--)
/*  looking at the integral between zeros[i] and zeros[i-1]
*/
    {
	mpfi_add(tmp,zeros[i],zeros[i-1]);
	mpfi_div_ui(tmp,tmp,2);
	mpfi_c_set_im(s,tmp);             // s <- 1/2 +(t1+t2)/2*I


	F(F_s,s);                        // F_s<-F(s)

/*
	printf("s is ");
	mpfi_c_print(s);
	printf("F(s) is ");
	mpfi_c_print(F_s);
*/    

	mpfi_sub(h,zeros[i],zeros[i-1]);
	mpfi_div_ui(h,h,2);

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
	    mpfi_c_div_z(c_tmp,ihk_n[n+n-1],nbang[n+n-1]);
	    mpfi_c_div_z(c_tmp1,ihk_n[n+n],nbang[n+n]);
	    mpfi_c_add(m_sum1_pos[n],m_sum1_pos[n-1],c_tmp1);
	    mpfi_c_add(m_sum1_neg[n],m_sum1_neg[n-1],c_tmp1);
	    mpfi_c_sub(m_sum1_pos[n],m_sum1_pos[n],c_tmp);
	    mpfi_c_add(m_sum1_neg[n],m_sum1_neg[n],c_tmp);
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
	    
	    mpfi_c_add_ui(c_tmp,ks,n+n+1);            // c_tmp<-ks+2n+1
	    mpfi_c_mul_z(c_tmp,c_tmp,nbang[n+n]);     // c_tmp<-(2n)!(ks+2n+1)
//	    printf("(2n)!(ks+2n+1) is ");mpfi_c_print(c_tmp);
	    mpfi_c_mul(c_tmp1,c_tmp,m_sum1_pos[n]);// c_tmp1<-(2n)!(ks+2n+1)M1+
	    mpfi_c_mul(c_tmp,c_tmp,m_sum1_neg[n]);  // c_tmp<-(2n)!(ks+2n+1)M1-
	    mpfi_c_sub(c_tmp1,c_tmp1,ihk_n[n+n+1]); // c_tmp1<-~-(ihk)^(2n+1)
	    mpfi_c_add(c_tmp,c_tmp,ihk_n[n+n+1]);   // c_tmp<-~+(ihk)^(2n+1)
	    mpfi_c_mul(c_tmp1,c_tmp1,m_sum2[n]);    
	    mpfi_c_mul(c_tmp,c_tmp,m_sum2[n]);
	    mpfi_c_div(c_tmp1,c_tmp1,ks_2n_kks[n]);
	    mpfi_c_div(c_tmp,c_tmp,ks_2n_kks[n]);
	    mpfi_c_add(n_sum_pos,n_sum_pos,c_tmp1);
	    mpfi_c_add(n_sum_neg,n_sum_neg,c_tmp);
//	    printf("iteration with n=%d\n",n);
//	    printf("m_sum1_pos is ");mpfi_c_print(m_sum1_pos[n]);
//	    printf("m_sum1_neg is ");mpfi_c_print(m_sum1_neg[n]);
//	    printf("m_sum2 is ");mpfi_c_print(m_sum2[n]);
//	    printf("pos adds ");mpfi_c_print(c_tmp1);
//	    printf("neg adds ");mpfi_c_print(c_tmp);
//	    printf("n_sum_pos is ");mpfi_c_print(n_sum_pos);
//	    printf("n_sum_neg is ");mpfi_c_print(n_sum_neg);
	};

	mpfi_c_exp(c_tmp,ihk_n[1]);
//	printf("exp(ihk) is ");mpfi_c_print(c_tmp);
	mpfi_c_mul(c_tmp,c_tmp,n_sum_pos);
//	printf("e(ihk)*n_sum_pos is ");mpfi_c_print(c_tmp);
	mpfi_c_neg(c_tmp1,ihk_n[1]);
	mpfi_c_exp(c_tmp1,c_tmp1);
//	printf("exp(-ihk) is ");mpfi_c_print(c_tmp1);
	mpfi_c_mul(c_tmp1,c_tmp1,n_sum_neg);
//	printf("e(-ihk)*n_sum_neg is ");mpfi_c_print(c_tmp1);
	mpfi_c_sub(c_tmp,c_tmp,c_tmp1);
	mpfi_c_mul(c_tmp,F_s,c_tmp);

//	printf("G(s+ih)-G(s-ih) at %d is ",i);
//	mpfi_c_print(c_tmp);

	mpfi_c_add(g_rho,g_rho,c_tmp);
	mpfi_c_conj(c_tmp,c_tmp);
	mpfi_c_add(g_rho,g_rho,c_tmp);


	mpfi_c_add(g_rho_sum,g_rho_sum,g_rho);

//	printf("g_rho of ");mpfi_c_print(s);
//	printf("      is ");mpfi_c_print(g_rho);
//	printf("sum   is ");mpfi_c_print(g_rho_sum);





	calc_taylor_error(g_rho_sum,s,h);




/*
	printf("zeros at height\n    ");
	mpfi_print(zeros[i]);
	printf("and ");
	mpfi_print(zeros[i-1]);
	printf("s is set to ");
	mpfi_c_print(s);
	printf("h is set to ");
	mpfi_print(h);
	printf("k*s is ");
	mpfi_c_print(ks);
	printf("i*h*k is ");
	mpfi_c_print(ihk);
*/	
	

    };
	
    mpfi_clear(half);
    mpfi_clear(h);
    mpfi_clear(tmp);
    mpfi_c_clear(s);
    mpfi_c_clear(ks);
    mpfi_c_clear(F_s);
    mpfi_c_clear(ks_sqr);
    mpfi_c_clear(n_sum_pos);
    mpfi_c_clear(n_sum_neg);
    mpfi_c_clear(c_tmp);
    mpfi_c_clear(c_tmp1);
    mpfi_c_clear(g_rho);
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
	printf("Usage: G <bits_prec> <x> <lambda> <taylor n>\n%lu <= bits_prec <= %lu\n",
	       MPFR_PREC_MIN,MPFR_PREC_MAX);
	printf("x >= 2\n");
	printf("lambda > 0\n");
	printf("taylor n > 0\n");
};

int main(int argc, char **argv)
{

    int prec;
    mpfi_t g1,ln2,err_term;
    mpfi_c_t g_rho_sum;

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


    if(argc!=5)
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
    
    mpfr_set_default_prec(prec);
    mpfi_c_setup(prec);


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


    printf("Lambda set to ");
    mpfr_print(lambda);
    printf("Taylor n set to %d.\n",taylor_n);


/* initialise our interval variables */
    printf("Initialising variables.\n");
    mpfi_c_init(g_rho_sum);
    mpfi_init(g1);
    mpfi_init(ln2);
    mpfi_init(err_term);


    mpfi_set_fr(lam_sqr,lambda);
    mpfi_mul(lam_sqr,lam_sqr,lam_sqr);
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

    printf("Entering g_sum.\n");
    if(g_sum(g_rho_sum)!=SUCCESS)
	return(QUIT_CODE);

    printf("Sum of G(rho) is ");mpfi_c_print(g_rho_sum);
/* not bothering to clear data structures on exit */    
    return(SUCCESS);
}
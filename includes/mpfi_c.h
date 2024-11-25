#ifndef MPFI_C
#define MPFI_C
//#include "mpfi.h"
//#include "mpfi_io.h"
#include "stdlib.h"
#include "math.h"
#define debug printf("Reached line %d.\n",__LINE__);
#define mpfi_inc(x,y) mpfi_add(x,x,y)
#define mpfi_inc_ui(x,y) mpfi_add_ui(x,x,y)
#define mpfi_dec(x,y) mpfi_sub(x,x,y)
#define mpfi_dec_ui(x,y) mpfi_sub_ui(x,x,y)
#define LNG_MIN_T (400000)
#define ZETA_ERR ((double) 1.02e-6) // assume t>10^6 and one term
#define FAILURE (1)

#ifndef __cplusplus
#define true (1==1)
#define false (1==0)
#else
extern "C" {
#endif

typedef struct{double re;double im;} complex;

inline complex cmplx (double re, double im)
{
  complex z;
  z.re=re;
  z.im=im;
  return(z);
}

inline complex c_neg (complex z)
{
  complex z1;
  z1.re=-z.re;
  z1.im=-z.im;
  return(z1);
}

void complex_print(complex s)
{
  printf("%10.8e+i%10.8e\n",s.re,s.im);
}

void mpfi_zero(mpfi_ptr res)
{
  mpfi_set_ui(res,0);
}


mpfi_t pow_res;
void mpfi_pow2 (mpfi_ptr x, mpfi_ptr y)
{
  mpfi_log(pow_res,x);
  mpfi_mul(pow_res,pow_res,y);
  mpfi_exp(x,pow_res);
}

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

void mpfi_c_print_str(const char *str, mpfi_c_ptr x)
{
  printf("%s",str);
  mpfi_c_print(x);
}

void mpfi_c_printn(mpfi_c_ptr z,int n)
{
    mpfi_out_str(stdout,10,n,z->re);
    printf(" + ");
    mpfi_out_str(stdout,10,n,z->im);
    printf("i\n");
};

void mpfi_print(mpfi_ptr x)
{
    mpfi_out_str(stdout,10,0,x);
    printf("\n");
};

void mpfi_print_str(const char *str, mpfi_ptr x)
{
  printf("%s",str);
  mpfi_print(x);
}

void mpfi_printn(mpfi_ptr x, int n)
{
    mpfi_out_str(stdout,10,n,x);
    printf("\n");
};


void mpfr_print(mpfr_ptr x)
{
    mpfr_out_str(stdout,10,0,x,GMP_RNDN);
    printf("\n");
};

void mpz_print(mpz_ptr x)
{
    mpz_out_str(stdout,10,x);
    printf("\n");
};

void mpfi_c_set(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_set(z1->re,z2->re);
    mpfi_set(z1->im,z2->im);
};

void mpfi_c_set_d(mpfi_c_ptr z, double re, double im)
{
  mpfi_set_d(z->re,re);
  mpfi_set_d(z->im,im);
};

void mpfi_c_set_ui(mpfi_c_ptr z, unsigned long int re, unsigned long int im)
{
  mpfi_set_ui(z->re,re);
  mpfi_set_ui(z->im,im);
};

void mpfi_c_set_i(mpfi_c_ptr z,mpfi_ptr x)
{
  mpfi_set(z->re,x);
  mpfi_set_ui(z->im,0);
}

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

void mpfi_c_add_i(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_ptr z3)
{
    mpfi_add(z1->re,z2->re,z3);
    mpfi_set(z1->im,z2->im);
};

#define mpfi_c_inc(x,y) mpfi_c_add(x,x,y)

void mpfi_c_sub(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_sub(z1->re,z2->re,z3->re);
    mpfi_sub(z1->im,z2->im,z3->im);
};

void mpfi_c_sub_i(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_ptr x)
{
    mpfi_sub(z1->re,z2->re,x);
    mpfi_set(z1->im,z2->im);
};

#define mpfi_c_dec(x,y) mpfi_c_sub(x,x,y)

void mpfi_c_dec_ui(mpfi_c_ptr z, unsigned long int i)
{
  mpfi_sub_ui(z->re,z->re,i);
};

void mpfi_c_sub_ui(mpfi_c_ptr res, mpfi_c_ptr z, unsigned long int i)
{
  mpfi_set(res->im,z->im);
  mpfi_sub_ui(res->re,z->re,i);
}

void mpfi_c_sub_d(mpfi_c_ptr res, mpfi_c_ptr z, double x)
{
  mpfi_set(res->im,z->im);
  mpfi_sub_d(res->re,z->re,x);
}

void mpfi_c_add_re(mpfi_c_ptr z1,mpfi_c_ptr z2,mpfi_ptr x)
{
    mpfi_add(z1->re,z2->re,x);
    mpfi_set(z1->im,z2->im);
};

void mpfi_c_add_ui(mpfi_c_ptr z1,mpfi_c_ptr z2, unsigned long int n)
{
    mpfi_add_ui(z1->re,z2->re,n);
    mpfi_set(z1->im,z2->im);
}

void mpfi_c_inc_ui(mpfi_c_ptr z,unsigned long int i)
{
  mpfi_add_ui(z->re,z->re,i);
}


void mpfi_c_neg(mpfi_c_ptr z1,mpfi_c_ptr z2)
{
    mpfi_neg(z1->im,z2->im);
    mpfi_neg(z1->re,z2->re);
};


mpfi_c_t d_tmp3,L_s;
mpfi_t m_tmp1,m_tmp2,m_tmp3,d_tmp1,d_tmp2,e_tmp,p_tmp1,p_tmp2;   
mpfi_t norm_tmp,c_log_tmp;
mpfr_t rel_err,rel_err1,erf_tmp1,erf_tmp2;
mpfi_t mod_tmp,mpfi_pi,mpfi_2_pi,mpfi_ln_2_pi,mpfi_sqrt_pi,mpfi_ln_pi;
mpfi_c_t lng_z,lng_tmp,lng_tmp1,lng_tmp2,ln_gamma_err,ln_gamma_err1;
mpfi_t log_2_pi_over_2;
mpfi_c_t c_one,c_minus_one,c_i,c_minus_i,c_zero;
mpfi_c_t lambda_tmp1;
mpfi_t lambda_tmp2;
mpfi_t mpfi_pi_by_2,mpfi_pi_by_4;
mpfi_t mpfi_log_2;

#define MAX_BERNOULLI_K (17)

mpfi_c_t hur_s_array[MAX_BERNOULLI_K],hur_tmp,hur_tmp1,hur_tmp2,hur_minus_z;
mpfi_t hur_mod,hur_alpha_n,pp_x_sigma,pp_tlnx,pow_tmp,hur_re_z;

mpfi_t bernoulli[MAX_BERNOULLI_K],h_bernoulli[MAX_BERNOULLI_K];

#define MIN_LNG_ABS (100)

mpfi_c_t theta_z;
mpfi_t zeta_term,zeta_temp,zeta_err,c0_temp,zeta_a;

void set_h_bernoulli()
{
  int i,j;
  for(i=0;i<MAX_BERNOULLI_K;i++)
    {
      mpfi_init(bernoulli[i]);
      mpfi_init(h_bernoulli[i]);
      mpfi_c_init(hur_s_array[i]);
    }
  mpfi_set_ui(bernoulli[0],1);
  mpfi_div_ui(bernoulli[0],bernoulli[0],6);
  mpfi_set_ui(bernoulli[1],1);
  mpfi_div_d(bernoulli[1],bernoulli[1],-30.0);		// b(4)
  mpfi_set_ui(bernoulli[2],1);
  mpfi_div_ui(bernoulli[2],bernoulli[2],42);	        // b(6)
  mpfi_set(bernoulli[3],bernoulli[1]);			// b(8)
  mpfi_set_ui(bernoulli[4],5);
  mpfi_div_ui(bernoulli[4],bernoulli[4],66);		// b(10)
  mpfi_set_d(bernoulli[5],-691.0);
  mpfi_div_ui(bernoulli[5],bernoulli[5],2730);		//b(12)
  mpfi_set_ui(bernoulli[6],7);
  mpfi_div_ui(bernoulli[6],bernoulli[6],6);		//b(14)
  mpfi_set_d(bernoulli[7],-3617);
  mpfi_div_ui(bernoulli[7],bernoulli[7],510);
  mpfi_set_ui(bernoulli[8],43867);
  mpfi_div_ui(bernoulli[8],bernoulli[8],798);
  mpfi_set_ui(bernoulli[9],174611);
  mpfi_div_ui(bernoulli[9],bernoulli[9],330);
  mpfi_neg(bernoulli[9],bernoulli[9]);
  mpfi_set_ui(bernoulli[10],854513);
  mpfi_div_ui(bernoulli[10],bernoulli[10],138);
  mpfi_set_d(bernoulli[11],-236364091.0);
  mpfi_div_ui(bernoulli[11],bernoulli[11],2730);
  mpfi_set_ui(bernoulli[12],8553103);
  mpfi_div_ui(bernoulli[12],bernoulli[12],6);
  mpfi_set_d(bernoulli[13],-23749461029.0);
  mpfi_div_ui(bernoulli[13],bernoulli[13],870);
  mpfi_set_ui(bernoulli[14],8615841276005);
  mpfi_div_ui(bernoulli[14],bernoulli[14],14322);
  mpfi_set_d(bernoulli[15],-7709321041217.0); // still < 2^53 so exact
  mpfi_div_ui(bernoulli[15],bernoulli[15],510);
  mpfi_set_ui(bernoulli[16],2577687858367);
  mpfi_div_ui(bernoulli[16],bernoulli[16],6);




  for(i=0;i<MAX_BERNOULLI_K;i++)
    {
      mpfi_set(h_bernoulli[i],bernoulli[i]);
      for(j=0;j<=i;j++)
	mpfi_div_ui(h_bernoulli[i],h_bernoulli[i],(j*2+1)*(j*2+2));
    }
  //mpfi_print(h_bernoulli[MAX_BERNOULLI_K-1]);

  //
  // Hare 1997 Prop 4.1 with \theta \in [0,\pi/2] so sec(\theta)\leq \sqrt{2}
  // max error is B_{2m}2^{2m}/2m/(2m-1)/|z|^{2m-1}
  mpfi_set(ln_gamma_err->re,bernoulli[MAX_BERNOULLI_K-1]);
  mpfi_mul_ui(ln_gamma_err->re,ln_gamma_err->re,1<<(MAX_BERNOULLI_K-1));
  mpfi_div_ui(ln_gamma_err->re,ln_gamma_err->re,MAX_BERNOULLI_K*(2*MAX_BERNOULLI_K-1));
  //mpfi_c_print(ln_gamma_err);
  mpfi_c_set(ln_gamma_err1,ln_gamma_err);
  for(i=1;i<2*MAX_BERNOULLI_K;i++)
    {
      mpfi_div_ui(ln_gamma_err->re,ln_gamma_err->re,MIN_LNG_ABS);
      mpfi_div_ui(ln_gamma_err1->re,ln_gamma_err1->re,LNG_MIN_T);
      //mpfi_c_print(ln_gamma_err);
    }
  mpfi_neg(ln_gamma_err->im,ln_gamma_err->re);
  //mpfi_c_print(ln_gamma_err);
  mpfi_put(ln_gamma_err->re,ln_gamma_err->im);
  //mpfi_c_print(ln_gamma_err);
  mpfi_set(ln_gamma_err->im,ln_gamma_err->re);

  mpfi_neg(ln_gamma_err1->im,ln_gamma_err1->re);
  //mpfi_c_print(ln_gamma_err1);
  mpfi_put(ln_gamma_err1->re,ln_gamma_err1->im);
  //mpfi_c_print(ln_gamma_err1);
  mpfi_set(ln_gamma_err1->im,ln_gamma_err1->re);
  //printf("Gamma Error set to ");mpfi_c_printn(ln_gamma_err1,10);

}

mpfr_t end_point,end_point1;


void mpfi_c_setup(int prec)
{
/* now initialise all global variables.
   we don't bother clearing them down afterwards, just exit. */
  int i;
    mpfr_set_default_prec(prec);
    mpfi_init(m_tmp1);
    mpfi_init(m_tmp2);
    mpfi_init(m_tmp3);
    mpfi_init(d_tmp1);
    mpfi_init(d_tmp2);
    mpfi_c_init(d_tmp3);
    mpfi_init(e_tmp);
    mpfi_init(p_tmp1);
    mpfi_init(p_tmp2);
    mpfi_init(pow_res);
    mpfi_init(norm_tmp);
    mpfr_init(rel_err);
    mpfr_init(rel_err1);
    mpfr_init(erf_tmp1);
    mpfr_init(erf_tmp2);
    mpfi_init(mod_tmp);
    mpfi_init(mpfi_pi);
    mpfi_const_pi(mpfi_pi);
    mpfi_init(mpfi_sqrt_pi);
    mpfi_sqrt(mpfi_sqrt_pi,mpfi_pi);
    mpfi_init(mpfi_2_pi);
    mpfi_init(mpfi_ln_2_pi);
    mpfi_mul_ui(mpfi_2_pi,mpfi_pi,2);
    mpfi_log(mpfi_ln_2_pi,mpfi_2_pi);
    mpfi_c_init(L_s);
    mpfi_c_init(lng_z);
    mpfi_c_init(lng_tmp);
    mpfi_c_init(lng_tmp1);
    mpfi_c_init(lng_tmp2);
    mpfi_c_init(ln_gamma_err);
    mpfi_c_init(ln_gamma_err1);
    mpfi_init(zeta_err);
    mpfi_set_d(zeta_err,ZETA_ERR);
    mpfi_neg(L_s->re,zeta_err);
    mpfi_put(zeta_err,L_s->re);
    mpfi_init(log_2_pi_over_2);
    mpfi_set(log_2_pi_over_2,mpfi_ln_2_pi);
    mpfi_div_ui(log_2_pi_over_2,log_2_pi_over_2,2);
    set_h_bernoulli();
    mpfi_init(c_log_tmp);
    mpfi_init(hur_mod);
    mpfi_init(hur_alpha_n);
    mpfi_c_init(hur_tmp);
    mpfi_c_init(hur_tmp1);
    mpfi_c_init(hur_tmp2);
    for(i=0;i<MAX_BERNOULLI_K;i++)
      mpfi_c_init(hur_s_array[i]);
    mpfi_c_init(hur_minus_z);
    mpfi_init(pp_tlnx);
    mpfi_init(pp_x_sigma);
    mpfi_init(pow_tmp);
    mpfi_c_init(c_one);
    mpfi_c_init(c_minus_one);
    mpfi_c_init(c_i);
    mpfi_c_init(c_minus_i);
    mpfi_c_init(c_zero);
    mpfi_c_set_ui(c_one,1,0);
    mpfi_c_neg(c_minus_one,c_one);
    mpfi_c_set_ui(c_i,0,1);
    mpfi_c_neg(c_minus_i,c_i);
    mpfi_c_set_ui(c_zero,0,0);
    mpfi_c_init(lambda_tmp1);
    mpfi_init(lambda_tmp2);
    mpfi_init(mpfi_pi_by_2);
    mpfi_div_ui(mpfi_pi_by_2,mpfi_pi,2);
    mpfi_init(mpfi_pi_by_4);
    mpfi_div_ui(mpfi_pi_by_4,mpfi_pi,4);
    mpfr_init(end_point);
    mpfr_init(end_point1);
    mpfi_init(mpfi_log_2);
    mpfi_const_log2(mpfi_log_2);
    mpfi_init(mpfi_ln_pi);
    mpfi_log(mpfi_ln_pi,mpfi_pi);
    mpfi_c_init(theta_z);
    mpfi_init(zeta_temp);
    mpfi_init(zeta_term);
    mpfi_init(c0_temp);
    mpfi_init(zeta_a);
    mpfi_init(hur_re_z);
}

void mpfi_c_norm(mpfi_ptr res, mpfi_c_ptr z)
{
  mpfi_sqr(res,z->re);
  mpfi_sqr(norm_tmp,z->im);
  mpfi_add(res,res,norm_tmp);
}

void mpfi_c_abs(mpfi_ptr res, mpfi_c_ptr z)
{
   mpfi_c_norm(res,z);
   mpfi_sqrt(res,res);
}


/* use the obvious method */

void mpfi_c_mul(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_mul(m_tmp1,z2->re,z3->re);
    mpfi_mul(m_tmp2,z2->im,z3->im);
    mpfi_sub(m_tmp3,m_tmp1,m_tmp2);
    mpfi_mul(m_tmp1,z2->re,z3->im);
    mpfi_mul(m_tmp2,z2->im,z3->re);
    mpfi_add(z1->im,m_tmp1,m_tmp2);
    mpfi_swap(z1->re,m_tmp3); // for cache reasons, might actually be quicker to use mpfi_set
};

void mpfi_c_sqr(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
  mpfi_sqr(m_tmp1,z2->re);
  mpfi_sqr(m_tmp2,z2->im);
  mpfi_mul(m_tmp3,z2->re,z2->im);
  mpfi_sub(z1->re,m_tmp1,m_tmp2);
  mpfi_mul_2ui(z1->im,m_tmp3,1);
}


/* use the M3 method */
/* This is actually slower and loses more precision ! */
/* are 2 adds and a sub slower than 1 mult? */
 /*
void mpfi_c_mul(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
  mpfi_add(m_tmp1,z2->re,z2->im); // (a+b)
  mpfi_add(m_tmp2,z3->re,z3->im); // (c+d)
  mpfi_mul(m_tmp3,m_tmp1,m_tmp2); // (a+b)(c+d)
  mpfi_mul(m_tmp1,z2->re,z3->re); // ac
  mpfi_mul(m_tmp2,z2->im,z3->im); // bd
  mpfi_sub(z1->re,m_tmp1,m_tmp2); // ac-bd
  mpfi_add(z1->im,m_tmp1,m_tmp2); // ac+bd
  mpfi_sub(z1->im,m_tmp3,z1->im); // (a+b)(c+d)-ac-bd
}
 */
/* multiply (mpfi_c_t z2) by (mpfi_t x) into (mpfi_c_t z1) safely */
void mpfi_c_mul_i(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_ptr x)
{
    mpfi_mul(z1->re,z2->re,x);
    mpfi_mul(z1->im,z2->im,x);
};

void mpfi_c_mul_d(mpfi_c_ptr z1, mpfi_c_ptr z2, double x)
{
    mpfi_mul_d(z1->re,z2->re,x);
    mpfi_mul_d(z1->im,z2->im,x);
}

void mpfi_c_mul_ui(mpfi_c_ptr z1, mpfi_c_ptr z2, unsigned long int x)
{
    mpfi_mul_ui(z1->re,z2->re,x);
    mpfi_mul_ui(z1->im,z2->im,x);
}

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
}

void mpfi_c_div_i (mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_ptr x)
{
	mpfi_div(z1->re,z1->re,x);
	mpfi_div(z1->im,z1->im,x);

}

void mpfi_c_div_ui (mpfi_c_ptr z1, mpfi_c_ptr z2, unsigned long int i)
{
    mpfi_div_ui(z1->re,z2->re,i);
    mpfi_div_ui(z1->im,z2->im,i);
};

void mpfi_c_div_d(mpfi_c_ptr z1, mpfi_c_ptr z2, double x)
{
  mpfi_div_d(z1->re,z2->re,x);
  mpfi_div_d(z1->im,z2->im,x);
}

void mpfi_c_exp(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_exp(e_tmp,z2->re);
    //debug;
//    mpfi_sin_cos(z1->im,z1->re,z2->im);   would be nice, but does not exist
    mpfi_cos(z1->re,z2->im);
    //debug;
    //mpfi_print_str("calling sin with ",z1->im);
    mpfi_sin(z1->im,z2->im);
    //debug;
    mpfi_c_mul_i(z1,z1,e_tmp);
};


void mpfi_c_pow_c(mpfi_c_ptr res, double i, complex s)
{
  //  printf("mpfi_c_pow_c called with %10.8e and ",i);
  //  mpfi_c_print(s);
  mpfi_set_d(p_tmp1,i);     
  mpfi_log(p_tmp1,p_tmp1);        // log(i)
  mpfi_mul_d(p_tmp2,p_tmp1,s.re);  // sigma*log(i)
  mpfi_exp(p_tmp2,p_tmp2);        // i^sigma
  mpfi_mul_d(p_tmp1,p_tmp1,s.im);  // t*log(i)
  mpfi_cos(res->re,p_tmp1);       // Re(res)=cos(t*log(i))
  //  printf("cos(t*log(i))=");
  //  mpfi_print(res->re);
  mpfi_sin(res->im,p_tmp1);       // Im(res)=sin(t*log(i))
  mpfi_c_mul_i(res,res,p_tmp2);   // res*=i^sigma
  //  printf("mpfi_c_pow returning ");
  //  mpfi_c_print(res);
}

void mpfi_c_pow_double_to_doubles(mpfi_c_ptr res, double x, double re_s, double im_s)
{
  mpfi_set_d(p_tmp1,x);     
  mpfi_log(p_tmp1,p_tmp1);        // log(x)
  mpfi_mul_d(p_tmp2,p_tmp1,re_s);  // sigma*log(x)
  mpfi_exp(p_tmp2,p_tmp2);        // x^sigma
  mpfi_mul_d(p_tmp1,p_tmp1,im_s);  // t*log(x)
  mpfi_cos(res->re,p_tmp1);       // Re(res)=cos(t*log(x))
  mpfi_sin(res->im,p_tmp1);       // Im(res)=sin(t*log(x))
  mpfi_c_mul_i(res,res,p_tmp2);   // res*=x^sigma
}


void mpfi_pow(mpfi_ptr res, mpfi_ptr x, mpfi_ptr y)
{
  mpfi_log(pow_tmp,x);
  mpfi_mul(pow_tmp,pow_tmp,y);
  mpfi_exp(res,pow_tmp);
}


void mpfi_c_pow_i_c (mpfi_c_ptr res, mpfi_ptr x, mpfi_c_ptr z)
{
  mpfi_pow(pp_x_sigma,x,z->re);
  //mpfi_print(pp_x_sigma);
  mpfi_log(pp_tlnx,x);
  //mpfi_print(pp_tlnx);
  mpfi_mul(pp_tlnx,pp_tlnx,z->im);
  //mpfi_print(pp_tlnx);
  mpfi_cos(res->re,pp_tlnx);
  mpfi_sin(res->im,pp_tlnx);
  mpfi_c_mul_i(res,res,pp_x_sigma);
}

/*
void mpfi_atan2 (mpfi_ptr res, mpfi_ptr y, mpfi_ptr x)
{
  mpfi_div(res,y,x);
  mpfi_atan(res,res);
  if(mpfi_is_neg(x))
    {
      if(mpfi_is_neg(y))
	mpfi_sub(res,res,mpfi_pi);
      else
	mpfi_add(res,res,mpfi_pi);
    }
}
*/

void mpfi_c_log (mpfi_c_ptr res, mpfi_c_ptr z)
{
  mpfi_c_norm(c_log_tmp,z);
  mpfi_log(c_log_tmp,c_log_tmp);
  mpfi_div_ui(c_log_tmp,c_log_tmp,2);
  //mpfi_div(res->im,z->im,z->re);
  //mpfi_print(res->im);
  mpfi_atan2(res->im,z->im,z->re);
  mpfi_set(res->re,c_log_tmp);
}

void mpfi_c_mod(mpfi_ptr res, mpfi_c_ptr z)
{
  mpfi_c_norm(res,z);
  mpfi_sqrt(res,res);
}

void mpfi_c_mod_c(mpfi_ptr res, complex z)
{
  mpfi_set_d(res,z.im);
  mpfi_sqr(res,res);
  mpfi_set_d(norm_tmp,z.re);
  mpfi_sqr(norm_tmp,norm_tmp);
  mpfi_add(res,res,norm_tmp);
  mpfi_sqrt(res,res);
}

int mpfi_rel_error(mpfi_ptr x)
{
  double err,dx;
  dx=fabs(mpfi_get_d(x));
  if(dx==0.0)
    return(999);
  mpfi_get_left(rel_err,x);
  mpfi_get_right(rel_err1,x);
  mpfr_sub(rel_err,rel_err,rel_err1,GMP_RNDN);
  err=fabs(mpfr_get_d(rel_err,GMP_RNDN));
  if(err==0.0)
    return(-999);
  return((int) log10(fabs(err/dx)));
}

int mpfi_abs_error(mpfi_ptr x)
{
  double err;
  mpfi_get_left(rel_err,x);
  mpfi_get_right(rel_err1,x);
  mpfr_sub(rel_err,rel_err,rel_err1,GMP_RNDN);
  err=fabs(mpfr_get_d(rel_err,GMP_RNDN));
  if(err==0.0)
    return(-999);
  return((int) log10(fabs(err)));
}


void mpfi_mod(mpfi_ptr res, mpfi_ptr den, mpfi_ptr num)
{
  int n;
  n=(int)floor(mpfi_get_d(den)/mpfi_get_d(num));
  mpfi_mul_ui(mod_tmp,num,n);
  mpfi_sub(res,den,mod_tmp);
}

void mpfi_c_arg(mpfi_ptr res, mpfi_c_ptr z)
{
  mpfi_div(res,z->im,z->re);
  mpfi_atan(res,res);
}

#define HI_BERNOULLI_K (9)

// quick and dirty version
// designed for z where Im(z) >> Re(z)
void mpfi_c_lngamma_hi(mpfi_c_ptr res, double re_z, double im_z)
{
  int n;
  if((re_z>=im_z)||(re_z<=0))
    {printf("mpfi_c_lngamma_hi called with re>=im and/or re<=0. Exiting.\n");
      exit(1);
    }
  //assert(re_z<im_z);
  //assert(re_z>0);
  mpfi_c_set_d(lng_z,re_z,im_z);
  mpfi_c_log(res,lng_z);
  mpfi_c_sub_d(lng_tmp,lng_z,0.5); // z-0.5
  //mpfi_c_print(lng_tmp);
  mpfi_c_mul(res,res,lng_tmp);    // (z-0.5)log(z)
  //mpfi_c_print(res);
  mpfi_c_sub(res,res,lng_z);      // (z-0.5)log(z)-z
  mpfi_c_add_i(res,res,log_2_pi_over_2); // (z-0.5)log(z)-z+log(2Pi)/2
  mpfi_c_set_i(lng_tmp,bernoulli[0]);  // B_2
  mpfi_c_div(lng_tmp,lng_tmp,lng_z);   // B_2/z
  mpfi_c_div_ui(lng_tmp,lng_tmp,2);    // B_2/2/z
  mpfi_c_add(res,res,lng_tmp);
  //mpfi_c_print(res);
  mpfi_c_mul(lng_tmp,lng_z,lng_z); //z^2
  n=2;
  while(n<=HI_BERNOULLI_K)
    {
      mpfi_c_mul(lng_z,lng_z,lng_tmp); //z^(2n-1)
      //mpfi_c_print(lng_z);
      mpfi_c_set_i(lng_tmp1,bernoulli[n-1]);
      mpfi_c_div(lng_tmp1,lng_tmp1,lng_z); // B_n/z^{2n-1}
      //mpfi_c_print(lng_tmp1);
      mpfi_c_div_ui(lng_tmp1,lng_tmp1,(n+n)*(n+n-1)); //B_n/z^{2n-1}/(2n(2n-1))
      //mpfi_c_print(lng_tmp1);
      mpfi_c_add(res,res,lng_tmp1);
      //mpfi_c_print(res);
      n++;
    }
  mpfi_div_d(lng_tmp->re,bernoulli[HI_BERNOULLI_K],2*re_z);
  mpfi_div_d(lng_tmp->re,lng_tmp->re,im_z);
  mpfi_div_ui(lng_tmp->re,lng_tmp->re,(HI_BERNOULLI_K*2+1)*(HI_BERNOULLI_K*2+2));
  for(n=1;n<=HI_BERNOULLI_K*2-1;n++)
    mpfi_div_d(lng_tmp->re,lng_tmp->re,im_z);
  mpfi_neg(lng_tmp->im,lng_tmp->re);
  mpfi_put(lng_tmp->re,lng_tmp->im);
  mpfi_set(lng_tmp->im,lng_tmp->re);
  

  //printf("lngamma error =");mpfi_c_printn(lng_tmp,20);
  mpfi_c_add(res,res,lng_tmp);
  return; 
}

void mpfi_c_lngamma(mpfi_c_ptr res, double re_z, double im_z)
{
  int n;
  if((re_z<0.0)||(im_z<0.0))
    {
      printf("lngamma only accepts arguments in first quadrant. Exiting.\n");
      exit(1);
    }
  if((re_z==0.0)&&(im_z==0.0))
    {
      printf("lngamma has a pole at 0. Exiting.\n");
      exit(0);
    }
  mpfi_c_set_d(lng_z,re_z,im_z);
  //mpfi_c_print(lng_z);
  if(im_z>=MIN_LNG_ABS)
  {
    mpfi_c_log(res,lng_z); // res=log(z)
    //mpfi_c_print(res);
    mpfi_c_sub_d(lng_tmp,lng_z,0.5); // z-0.5
    //mpfi_c_print(lng_tmp);
    mpfi_c_mul(res,res,lng_tmp);    // (z-0.5)log(z)
    //mpfi_c_print(res);
    mpfi_c_sub(res,res,lng_z);      // (z-0.5)log(z)-z
    mpfi_c_add_i(res,res,log_2_pi_over_2); // (z-0.5)log(z)-z+log(2Pi)/2
    mpfi_c_set_i(lng_tmp,bernoulli[0]);  // B_2
    mpfi_c_div(lng_tmp,lng_tmp,lng_z);   // B_2/z
    mpfi_c_div_ui(lng_tmp,lng_tmp,2);    // B_2/2/z
    mpfi_c_add(res,res,lng_tmp);
    //mpfi_c_print(res);
    mpfi_c_mul(lng_tmp,lng_z,lng_z); //z^2
    n=2;
    while(n<=MAX_BERNOULLI_K)
      {
	mpfi_c_mul(lng_z,lng_z,lng_tmp); //z^(2n-1)
	//mpfi_c_print(lng_z);
	mpfi_c_set_i(lng_tmp1,bernoulli[n-1]);
	mpfi_c_div(lng_tmp1,lng_tmp1,lng_z); // B_n/z^{2n-1}
	//mpfi_c_print(lng_tmp1);
	mpfi_c_div_ui(lng_tmp1,lng_tmp1,(n+n)*(n+n-1)); //B_n/z^{2n-1}/(2n(2n-1))
	//mpfi_c_print(lng_tmp1);
	mpfi_c_add(res,res,lng_tmp1);
	//mpfi_c_print(res);
	n++;
      }
    //mpfi_c_print(ln_gamma_err);
    mpfi_c_add(res,res,ln_gamma_err);
    return;
  }
  else
    {
      
      mpfi_c_log(lng_tmp2,lng_z);
      for(n=1;n<MIN_LNG_ABS;n++)
	{
	  mpfi_c_add_ui(lng_z,lng_z,1);
	  mpfi_c_log(lng_tmp1,lng_z);
	  mpfi_c_add(lng_tmp2,lng_tmp2,lng_tmp1);
	}
      mpfi_c_add_ui(lng_z,lng_z,1); // z+10
	  mpfi_c_log(res,lng_z); // res=log(z)
	  //mpfi_c_print(res);
	  mpfi_c_sub_d(lng_tmp,lng_z,0.5); // z-0.5
	  //mpfi_c_print(lng_tmp);
	  mpfi_c_mul(res,res,lng_tmp);    // (z-0.5)log(z)
	  //mpfi_c_print(res);
	  mpfi_c_sub(res,res,lng_z);      // (z-0.5)log(z)-z
	  mpfi_c_add_i(res,res,log_2_pi_over_2); // (z-0.5)log(z)-z+log(2Pi)/2
	  mpfi_c_set_i(lng_tmp,bernoulli[0]);  // B_2
	  mpfi_c_div(lng_tmp,lng_tmp,lng_z);   // B_2/z
	  mpfi_c_div_ui(lng_tmp,lng_tmp,2);    // B_2/2/z
	  mpfi_c_add(res,res,lng_tmp);
	  //mpfi_c_print(res);
	  mpfi_c_mul(lng_tmp,lng_z,lng_z); //z^2
	  n=2;
	  while(n<=MAX_BERNOULLI_K)
	  {
		  mpfi_c_mul(lng_z,lng_z,lng_tmp); //z^(2n-1)
		  //mpfi_c_print(lng_z);
		  mpfi_c_set_i(lng_tmp1,bernoulli[n-1]);
		  mpfi_c_div(lng_tmp1,lng_tmp1,lng_z); // B_n/z^{2n-1}
		  //mpfi_c_print(lng_tmp1);
		  mpfi_c_div_ui(lng_tmp1,lng_tmp1,(n+n)*(n+n-1)); //B_n/z^{2n-1}/(2n(2n-1))
		  //mpfi_c_print(lng_tmp1);
		  mpfi_c_add(res,res,lng_tmp1);
		  //mpfi_c_print(res);
		  n++;
	  }
	  mpfi_c_sub(res,res,lng_tmp2);
	  //mpfi_c_print(ln_gamma_err);
	  mpfi_c_add(res,res,ln_gamma_err);
	  return;
  }
}

void mpfi_c_gamma(mpfi_c_ptr res, double re_z, double im_z)
{
  mpfi_c_lngamma(res,re_z,im_z);
  mpfi_c_exp(res,res);
}

#define DEFAULT_HUR_N (50)

//
// set hur_s_array[j] to B_{2j+2}s(s+1)..(s+2j)/(2j+2)!
// for j=0..MAX_BERNOULLI_K-1
//
void mpfi_c_hur_init1();

void mpfi_c_hur_init(double re_z, double im_z)
{
	mpfi_c_set_d(hur_tmp,re_z,im_z);
	mpfi_c_hur_init1();
}

void mpfi_c_hur_init1()
{
  long int i;
	mpfi_c_set(hur_s_array[0],hur_tmp); // = z
	for(i=1;i<MAX_BERNOULLI_K;i++)
	{
		mpfi_c_add_ui(hur_tmp,hur_tmp,1);  // z+i
		mpfi_c_mul(hur_s_array[i],hur_s_array[i-1],hur_tmp);
		mpfi_c_add_ui(hur_tmp,hur_tmp,1);  // z+i+1
		mpfi_c_mul(hur_s_array[i],hur_s_array[i],hur_tmp); 
	}
	for(i=0;i<MAX_BERNOULLI_K;i++)
	{
		//mpfi_print(h_bernoulli[i]);
		//mpfi_c_print(hur_s_array[i]);
		mpfi_c_mul_i(hur_s_array[i],hur_s_array[i],h_bernoulli[i]);
		//printf("hur_s_array[%d]=",i);mpfi_c_print(hur_s_array[i]);
	}
}

void mpfi_c_hurwitz2(mpfi_c_ptr, mpfi_ptr, long int);

int N_factor=2;
void mpfi_c_hurwitz1(mpfi_c_ptr res, double re_z, double im_z, mpfi_ptr alpha)
{
  int N;
  N=ceil(im_z)*N_factor; // this will do for modulus
  if(N<DEFAULT_HUR_N)
    N=DEFAULT_HUR_N; // n is how many terms to sum
  mpfi_c_set_d(hur_minus_z,-re_z,-im_z); // -z
  mpfi_c_hurwitz2(res,alpha,N);
}

// hur_minus_z = -z on entry
void mpfi_c_hurwitz2(mpfi_c_ptr res,mpfi_ptr alpha, long int N)
{
  long int i;
  double re_z;

  mpfi_set(hur_alpha_n,alpha);           
  mpfi_c_pow_i_c(res,hur_alpha_n,hur_minus_z);   // res = alpha^(-z) 
  //mpfi_c_print(res);
  for(i=1;i<N;i++)
    {
      mpfi_add_ui(hur_alpha_n,hur_alpha_n,1);      
      mpfi_c_pow_i_c(hur_tmp,hur_alpha_n,hur_minus_z); // (alpha+i)^(-z)
      mpfi_c_add(res,res,hur_tmp);
    }
  //printf("alpha=");mpfi_print(alpha);
  //printf("sum n=0..%d (n+alpha)^(-s)=",N-1);mpfi_c_print(res);
  //exit(0);
  // res = sum_{m=0..N-1} (m+alpha)^(-s)
  mpfi_add_ui(hur_alpha_n,alpha,N); // hur_alpha_n=(N+alpha)
  mpfi_c_add_ui(hur_tmp1,hur_minus_z,1); // hur_tmp1=-s+1
  mpfi_c_pow_i_c(hur_tmp,hur_alpha_n,hur_tmp1); // hur_tmp=(N+alpha)^(-s+1)
  mpfi_c_div(hur_tmp2,hur_tmp,hur_tmp1);    // hur_tmp2=-(N+alpha)^(-s+1)/(s-1)
  mpfi_c_sub(res,res,hur_tmp2);             // res=sum + (N+alpha)^(-s+1)/(s-1)
  //printf("sum+(%d+alpha)^(-s+1)/(s-1)=",N);mpfi_c_print(res); 
  //exit(0);
  mpfi_c_div_i(hur_tmp,hur_tmp,hur_alpha_n); // hur_tmp=(N+alpha)^(-s)
  mpfi_c_div_ui(hur_tmp2,hur_tmp,2);         // hur_tmp2=(N+alpha)^(-s)/2
  mpfi_c_add(res,res,hur_tmp2);              // res=sum + (N+alpha)^(-s+1)/(s-1) + (N+alpha)^(-s)/2
  //printf("sum+(%d+alpha)^(-s+1)/(s-1)+(%d+alpha)^(-s)/2=",N,N);mpfi_c_print(res);
  //exit(0);
  mpfi_c_div_i(hur_tmp,hur_tmp,hur_alpha_n); // (N+alpha)^(-s-1)
  mpfi_c_mul(hur_tmp1,hur_s_array[0],hur_tmp); // B_2s(N+alpha)^(-s-1)/2!
  mpfi_c_add(res,res,hur_tmp1);
  //printf("sum+(%d+alpha)^(-s+1)/(s-1)+(%d+alpha)^(-s)/2+B_2s(%d+alpha)^(-s-1)/2!=\n",N,N,N);mpfi_c_print(res);
  //exit(0);
  mpfi_mul(hur_alpha_n,hur_alpha_n,hur_alpha_n); // (N+alpha)^2
  //printf("(N+alpha)^2=");mpfi_print(hur_alpha_n);
  for(i=1;i<MAX_BERNOULLI_K-1;i++)
    {
		mpfi_c_div_i(hur_tmp,hur_tmp,hur_alpha_n); // (N+alpha)^(-s-2i-1)
		mpfi_c_mul(hur_tmp1,hur_tmp,hur_s_array[i]);
		mpfi_c_inc(res,hur_tmp1);
    }
  mpfi_c_div_i(hur_tmp,hur_tmp,hur_alpha_n); // (N+alpha)^(-s-2i-1)
  mpfi_c_mul(hur_tmp1,hur_tmp,hur_s_array[i]);
  mpfi_c_inc(res,hur_tmp1);
  
  // hur_tmp1 = first term omitted
  mpfi_c_neg(hur_minus_z,hur_minus_z); // now = z
  mpfi_c_add_ui(hur_minus_z,hur_minus_z,(2*MAX_BERNOULLI_K-1));
  mpfi_c_mod(hur_mod,hur_minus_z);
  mpfi_c_mod(hur_alpha_n,hur_tmp1);
  mpfi_mul(hur_mod,hur_mod,hur_alpha_n);
  mpfi_sub_ui(hur_re_z,hur_minus_z->re,2*MAX_BERNOULLI_K-1);
  mpfi_neg(hur_re_z,hur_re_z);
  mpfi_div(hur_alpha_n,hur_mod,hur_re_z);//+2*MAX_BERNOULLI_K-1);
  //printf("hur=");mpfi_c_print(res);
  //printf("error=");mpfi_print(hur_alpha_n);
  mpfi_neg(hur_mod,hur_alpha_n);
  mpfi_put(hur_alpha_n,hur_mod);
  //printf("error=");mpfi_print(hur_alpha_n);
  mpfi_inc(res->re,hur_alpha_n);
  mpfi_inc(res->im,hur_alpha_n);
  //printf("hur=");mpfi_c_print(res);
  //exit(0);
}

void mpfi_c_hurwitz(mpfi_c_ptr res, double re_z, double im_z, mpfi_ptr alpha)
{
  if((re_z<0.0)||(im_z<0.0)||(!mpfi_is_strictly_pos(alpha)))
    {
      printf("hurwitz only works with z in 1st quadrant and alpha>0.0. Exiting.\n");
      exit(1);
    }
  mpfi_c_hur_init(re_z,im_z);
  mpfi_c_hurwitz1(res,re_z,im_z,alpha);
}

void mpfi_c_hurwitz_c(mpfi_c_ptr res, mpfi_c_ptr z, mpfi_ptr alpha)
{
  long int N;
  if((mpfi_cmp_d(z->re,0.0)<=0)||(mpfi_cmp_d(z->im,0.0)<=0)||(mpfi_cmp_d(alpha,0.0)<=0))
    {
      printf("hurwitz only works with z in 1st quadrant and alpha>0.0. Exiting.\n");
      exit(1);
    }
  mpfi_c_set(hur_tmp,z);
  mpfi_c_hur_init1();
  mpfi_c_neg(hur_minus_z,z);
  N=mpfi_get_d(z->im)*N_factor;
  if(N<DEFAULT_HUR_N)
    N=DEFAULT_HUR_N; // N is how many terms to sum
  mpfi_c_hurwitz2(res,alpha,N);
}

void write_mpfi (FILE *outfile, mpfi_ptr z)
{
  double x[2];
  mpfi_get_left(end_point,z);
  x[0]=mpfr_get_d(end_point,GMP_RNDD);
  mpfi_get_right(end_point,z);
  x[1]=mpfr_get_d(end_point,GMP_RNDU);
  if(fwrite(x,sizeof(double),2,outfile)!=2)
    {
      printf("Error writing data in write_mpfi. Exiting.\n");
      exit(0);
    }
}

// convert complex interval to 4 doubles and output.
void write_mpfi_c (FILE *outfile, mpfi_c_ptr z)
{
    double x[4];

    mpfi_get_left(end_point,z->re);
    x[0]=mpfr_get_d(end_point,GMP_RNDD);
    mpfi_get_right(end_point,z->re);
    x[1]=mpfr_get_d(end_point,GMP_RNDU);
    mpfi_get_left(end_point,z->im);
    x[2]=mpfr_get_d(end_point,GMP_RNDD);
    mpfi_get_right(end_point,z->im);
    x[3]=mpfr_get_d(end_point,GMP_RNDU);
    if(fwrite(x,sizeof(double),4,outfile)!=4)
      {
	printf("Error writing data in write_mpfi_c. Exiting.\n");
	exit(0);
      }
}

void mpfr_write_bin(FILE *outfile, mpfr_ptr x)
{
  fwrite(x,sizeof(mpfr_t),1,outfile);
  int num_limbs=ceil(x->_mpfr_prec/mp_bits_per_limb);
  fwrite(x->_mpfr_d,sizeof(mp_limb_t),num_limbs,outfile);
}

void mpfi_write_bin(FILE *outfile, mpfi_ptr z)
{
  mpfi_get_left(end_point,z);
  mpfr_write_bin(outfile,end_point);
  mpfi_get_right(end_point,z);
  mpfr_write_bin(outfile,end_point);
}

void mpfr_read_bin(FILE *infile, mpfr_ptr x)
{
  mp_prec_t prec=mpfr_get_prec(x);
  fread(x,sizeof(mpfr_t),1,infile);
  if(x->_mpfr_prec!=prec)
    {printf("Precision mismatch in mpfr_read_bin. Exiting.\n");exit(0);}
  int num_limbs=ceil(x->_mpfr_prec/mp_bits_per_limb);
  fread(x->_mpfr_d,sizeof(mp_limb_t),num_limbs,infile);
}


void mpfi_read_bin(FILE *infile, mpfi_ptr z)
{
  mpfr_read_bin(infile,end_point);
  mpfr_read_bin(infile,end_point1);
  mpfi_interv_fr(z,end_point,end_point1);
}

int mpfi_contains_zero(mpfi_ptr x)
{
  return(mpfi_cmp_d(x,0.0)==0);
}

int mpfi_c_contains_zero(mpfi_c_ptr z)
{
  return(mpfi_contains_zero(z->re)&&mpfi_contains_zero(z->im));
}

void mpfi_c_zero(mpfi_c_ptr z)
{
   mpfi_set_ui(z->re,0);
   mpfi_set_ui(z->im,0);
}

void mpfi_c_one(mpfi_c_ptr z)
{
  mpfi_set_ui(z->re,1);
  mpfi_set_ui(z->im,0);
}

void mpfi_c_sqrt(mpfi_c_ptr res, mpfi_c_ptr z)
{
   mpfi_c_log(res,z);
   mpfi_c_div_ui(res,res,2);
   mpfi_c_exp(res,res);
   if(mpfi_is_neg(res->re))
     mpfi_c_neg(res,res);
}

void mpfi_erf(mpfi_ptr res, mpfi_ptr op1)

{
// a nice strictly increasing function so nothing clever required
      mpfi_get_left(erf_tmp1,op1);
      mpfr_erf(erf_tmp1,erf_tmp1,GMP_RNDD);
      mpfi_get_right(erf_tmp2,op1);
      mpfr_erf(erf_tmp2,erf_tmp2,GMP_RNDU);
      mpfi_interv_fr(res,erf_tmp1,erf_tmp2);
};

void mpfi_erfc(mpfi_ptr res, mpfi_ptr x)
{
  mpfi_erf(res,x);
  mpfi_neg(res,res);
  mpfi_add_ui(res,res,1);
}

// Im log Gamma(0.25+it/2)-tlog(Pi)/2
void mpfi_theta(mpfi_t res, double t)
{
  double t_2=t/2.0;
  mpfi_c_lngamma_hi(theta_z,0.25,t_2);
  mpfi_mul_d(res,mpfi_ln_pi,-t_2);
  mpfi_add(res,res,theta_z->im);
}

void mpfi_c0(mpfi_ptr res, mpfi_ptr rho)
{
  mpfi_mul(c0_temp,rho,mpfi_2_pi);
  mpfi_cos(c0_temp,c0_temp);
  mpfi_sqr(res,rho);
  mpfi_sub(res,res,rho);
  mpfi_sub_d(res,res,1.0/16.0);
  mpfi_mul(res,res,mpfi_2_pi);
  mpfi_cos(res,res);
  mpfi_div(res,res,c0_temp);
}

// calculate zeta(1/2+it) using Riemann Siegel
// with 10 terms
// error = 25966*t^(-23/4)
// logs[i]=log(i+1) i=0..sqrt(t/2Pi)
// sqrts[i]=1/sqrt(i+1)
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void mpfi_zeta_rs(mpfi_ptr res, double t, mpfi_t *logs, mpfi_t *sqrts)
{
  int nu,k;

  if(t<1000000.0)
    {
      printf("Fatal error in mpfi_c_zeta_rs. t too small (<10^6). Exiting.\n");
      exit(0);
    }

  nu=sqrt(t/2.0/M_PI);
  mpfi_zero(res);
  mpfi_theta(zeta_temp,t);
  for(k=0;k<nu;k++)
    {
      mpfi_mul_d(zeta_term,logs[k],-t);
      mpfi_add(zeta_term,zeta_term,zeta_temp);
      mpfi_cos(zeta_term,zeta_term);
      mpfi_mul(zeta_term,zeta_term,sqrts[k]);
      mpfi_add(res,res,zeta_term);
    }
  mpfi_mul_ui(res,res,2); // the main term
  mpfi_set_d(zeta_a,t/2.0);
  mpfi_div(zeta_a,zeta_a,mpfi_pi);
  mpfi_sqrt(zeta_a,zeta_a);         // a
  mpfi_sub_ui(zeta_temp,zeta_a,nu); // temp=rho
  mpfi_c0(zeta_term,zeta_temp); // term=c0(rho)
  //printf("c0(rho)=");mpfi_printn(zeta_term,30);
  mpfi_sqrt(zeta_temp,zeta_a);
  mpfi_div(zeta_term,zeta_term,zeta_temp); // c0(rho)/sqrt(a)
  if(nu&1)
    mpfi_add(res,res,zeta_term);
  else
    mpfi_sub(res,res,zeta_term);
  mpfi_add(res,res,zeta_err);
}

// convert a 128 bit unsigned int into mpfi_t
void mpfi_set_128(mpfi_ptr res, __uint128_t x)
{
  __uint64_t tmp=(__uint64_t) x; // just gets bottom 64 bits
  mpfi_set_ui(res,x>>64); // gets top 64 bits
  mpfi_mul_2ui(res,res,64);
  mpfi_add_ui(res,res,tmp);
}


/*

Stuff below here got superceded

// This was only used by windowed/g_even.c which
// I think is redundant
//
void mpfi_c_lngc(mpfi_c_ptr res, mpfi_c_ptr z)
{
  int n;
  mpfi_c_set(lng_z,z);
  if((!mpfi_is_strictly_pos(lng_z->re))||(!mpfi_is_strictly_pos(lng_z->im)))
    {
      printf("lngamma only accepts arguments in first quadrant. Exiting.\n");
      exit(1);
    }
  //mpfi_c_print(lng_z);
  if(mpfi_cmp_d(lng_z->im,LNG_MIN_T)>0)
    {
      mpfi_c_log(res,lng_z); // res=log(z)
      //mpfi_c_print(res);
      mpfi_c_sub_d(lng_tmp,lng_z,0.5); // z-0.5
      //mpfi_c_print(lng_tmp);
      mpfi_c_mul(res,res,lng_tmp);    // (z-0.5)log(z)
      //mpfi_c_print(res);
      mpfi_c_sub(res,res,lng_z);      // (z-0.5)log(z)-z
      mpfi_c_add_i(res,res,log_2_pi_over_2); // (z-0.5)log(z)-z+log(2Pi)/2
      mpfi_c_set_i(lng_tmp,bernoulli[0]);  // B_2
      mpfi_c_div(lng_tmp,lng_tmp,lng_z);   // B_2/z
      mpfi_c_div_ui(lng_tmp,lng_tmp,2);    // B_2/2/z
      mpfi_c_add(res,res,lng_tmp);
      //mpfi_c_print(res);
      mpfi_c_mul(lng_tmp,lng_z,lng_z); //z^2
      n=2;
      while(n<=MAX_BERNOULLI_K)
	{
	  mpfi_c_mul(lng_z,lng_z,lng_tmp); //z^(2n-1)
	  //mpfi_c_print(lng_z);
	  mpfi_c_set_i(lng_tmp1,bernoulli[n-1]);
	  mpfi_c_div(lng_tmp1,lng_tmp1,lng_z); // B_n/z^{2n-1}
	  //mpfi_c_print(lng_tmp1);
	  mpfi_c_div_ui(lng_tmp1,lng_tmp1,(n+n)*(n+n-1)); //B_n/z^{2n-1}/(2n(2n-1))
	  //mpfi_c_print(lng_tmp1);
	  mpfi_c_add(res,res,lng_tmp1);
	  //mpfi_c_print(res);
	  n++;
	}
      //printf("lng error=");mpfi_c_print(ln_gamma_err1);
      mpfi_c_add(res,res,ln_gamma_err1);
      return;
    }
  else
    {
      printf("This version of lnGamma requires Im(z)>%d. Exiting.\n",LNG_MIN_T);mpfi_c_print(z);
      exit(FAILURE); 
      mpfi_c_log(lng_tmp2,lng_z);
      for(n=1;n<10;n++)
	{
	  mpfi_c_add_ui(lng_z,lng_z,1);
	  mpfi_c_log(lng_tmp1,lng_z);
	  mpfi_c_add(lng_tmp2,lng_tmp2,lng_tmp1);
	}
      mpfi_c_add_ui(lng_z,lng_z,1); // z+10
      mpfi_c_log(res,lng_z); // res=log(z)
      //mpfi_c_print(res);
      mpfi_c_sub_d(lng_tmp,lng_z,0.5); // z-0.5
      //mpfi_c_print(lng_tmp);
      mpfi_c_mul(res,res,lng_tmp);    // (z-0.5)log(z)
      //mpfi_c_print(res);
      mpfi_c_sub(res,res,lng_z);      // (z-0.5)log(z)-z
      mpfi_c_add_i(res,res,log_2_pi_over_2); // (z-0.5)log(z)-z+log(2Pi)/2
      mpfi_c_set_i(lng_tmp,bernoulli[0]);  // B_2
      mpfi_c_div(lng_tmp,lng_tmp,lng_z);   // B_2/z
      mpfi_c_div_ui(lng_tmp,lng_tmp,2);    // B_2/2/z
      mpfi_c_add(res,res,lng_tmp);
      //mpfi_c_print(res);
      mpfi_c_mul(lng_tmp,lng_z,lng_z); //z^2
      n=2;
      while(n<=MAX_BERNOULLI_K)
	{
	  mpfi_c_mul(lng_z,lng_z,lng_tmp); //z^(2n-1)
	  //mpfi_c_print(lng_z);
	  mpfi_c_set_i(lng_tmp1,bernoulli[n-1]);
	  mpfi_c_div(lng_tmp1,lng_tmp1,lng_z); // B_n/z^{2n-1}
	  //mpfi_c_print(lng_tmp1);
	  mpfi_c_div_ui(lng_tmp1,lng_tmp1,(n+n)*(n+n-1)); //B_n/z^{2n-1}/(2n(2n-1))
	  //mpfi_c_print(lng_tmp1);
	  mpfi_c_add(res,res,lng_tmp1);
	  //mpfi_c_print(res);
	  n++;
	}
      mpfi_c_sub(res,res,lng_tmp2);
      //mpfi_c_print(ln_gamma_err);
      mpfi_c_add(res,res,ln_gamma_err);
      return;
    }
}
*/
#ifdef __cplusplus
}
#endif
#endif

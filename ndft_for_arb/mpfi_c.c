#include "stdio.h"
#include "mpfi.h"
#include "mpfi_c.h"

// First the mpfi_c stuff

/* initialisation */
inline void mpfi_c_init(mpfi_c_ptr z)
{
    mpfi_init(z->re);
    mpfi_init(z->im);
};

/* clearing */
inline void mpfi_c_clear(mpfi_c_ptr z)
{
    mpfi_clear(z->re);
    mpfi_clear(z->im);
};

/* efficient swapping */
inline void mpfi_c_swap(mpfi_c_ptr z1,mpfi_c_ptr z2)
{
    mpfi_swap(z1->re,z2->re);
    mpfi_swap(z1->im,z2->im);
};

/* simple printing */
inline void mpfi_c_print(mpfi_c_ptr z)
{
    mpfi_out_str(stdout,10,0,z->re);
    printf(" + ");
    mpfi_out_str(stdout,10,0,z->im);
    printf("i\n");
};

inline void mpfi_c_print_str(const char *str, mpfi_c_ptr x)
{
  printf("%s",str);
  mpfi_c_print(x);
}

inline void mpfi_c_printn(mpfi_c_ptr z,int n)
{
    mpfi_out_str(stdout,10,n,z->re);
    printf(" + ");
    mpfi_out_str(stdout,10,n,z->im);
    printf("i\n");
};

inline void mpfi_print(mpfi_ptr x)
{
    mpfi_out_str(stdout,10,0,x);
    printf("\n");
};

inline void mpfi_print_str(const char *str, mpfi_ptr x)
{
  printf("%s",str);
  mpfi_print(x);
}

inline void mpfi_printn(mpfi_ptr x, int n)
{
    mpfi_out_str(stdout,10,n,x);
    printf("\n");
};

// setting from various data types
inline void mpfi_c_set(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_set(z1->re,z2->re);
    mpfi_set(z1->im,z2->im);
};

inline void mpfi_c_set_d(mpfi_c_ptr z, double re, double im)
{
  mpfi_set_d(z->re,re);
  mpfi_set_d(z->im,im);
};

inline void mpfi_c_set_ui(mpfi_c_ptr z, unsigned long int re, unsigned long int im)
{
  mpfi_set_ui(z->re,re);
  mpfi_set_ui(z->im,im);
};

// just set the real part
inline void mpfi_c_set_re(mpfi_c_ptr z,mpfi_ptr x)
{
    mpfi_set(z->re,x);
};

inline void mpfi_c_set_re_d(mpfi_c_ptr z,double x)
{
    mpfi_set_d(z->re,x);
};

// just set the imaginary part
inline void mpfi_c_set_im(mpfi_c_ptr z, mpfi_ptr x)
{
    mpfi_set(z->im,x);
};

// z1=z2+z3
inline void mpfi_c_add(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_add(z1->re,z2->re,z3->re);
    mpfi_add(z1->im,z2->im,z3->im);
};

// z1=z2-z3
inline void mpfi_c_sub(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_sub(z1->re,z2->re,z3->re);
    mpfi_sub(z1->im,z2->im,z3->im);
};

// z1=-z2
inline void mpfi_c_neg(mpfi_c_ptr z1,mpfi_c_ptr z2)
{
    mpfi_neg(z1->im,z2->im);
    mpfi_neg(z1->re,z2->re);
};

// globally defined variables
// used to avoid repeated init/clear
//
mpfi_c_t d_tmp3,d_tmp5,c_sqrt_tmp;
mpfi_t m_tmp1,m_tmp2,m_tmp3,d_tmp1,d_tmp2,d_tmp4,e_tmp,p_tmp1,p_tmp2;   
mpfi_t norm_tmp,norm_tmp1,c_log_tmp,c_log_tmp1,c_arg_tmp;
mpfi_t mod_tmp,mpfi_ln_2_pi,mpfi_sqrt_pi,mpfi_ln_pi;
mpfi_c_t lng_z,lng_tmp,lng_tmp1,lng_tmp2,ln_gamma_err,ln_gamma_err1;
mpfi_t log_2_pi_over_2;
mpfi_c_t c_one,c_minus_one,c_i,c_minus_i,c_zero;

mpfi_t mpfi_pi_by_2,mpfi_pi_by_4;
mpfi_t mpfi_log_2;


mpfi_t pp_x_sigma,pp_tlnx,pow_tmp,pow_tmp1;


inline void mpfi_c_setup(unsigned long int prec)
{
/* now initialise all global variables.
   we don't bother clearing them down afterwards, just exit. */
  mpfr_set_default_prec(prec);
  mpfi_init(m_tmp1);
  mpfi_init(m_tmp2);
  mpfi_init(m_tmp3);
  mpfi_init(d_tmp1);
  mpfi_init(d_tmp2);
  mpfi_init(d_tmp4);
  mpfi_c_init(d_tmp3);
  mpfi_c_init(d_tmp5);
  mpfi_init(e_tmp);
  mpfi_init(p_tmp1);
  mpfi_init(p_tmp2);
  mpfi_init(norm_tmp);
  mpfi_init(norm_tmp1);
  mpfi_init(mod_tmp);
  mpfi_init(mpfi_pi);
  mpfi_const_pi(mpfi_pi);
  mpfi_init(mpfi_sqrt_pi);
  mpfi_sqrt(mpfi_sqrt_pi,mpfi_pi);
  mpfi_init(mpfi_2_pi);
  mpfi_init(mpfi_ln_2_pi);
  mpfi_mul_ui(mpfi_2_pi,mpfi_pi,2);
  mpfi_log(mpfi_ln_2_pi,mpfi_2_pi);
  mpfi_init(log_2_pi_over_2);
  mpfi_set(log_2_pi_over_2,mpfi_ln_2_pi);
  mpfi_div_ui(log_2_pi_over_2,log_2_pi_over_2,2);
  mpfi_init(c_log_tmp);
  mpfi_init(c_log_tmp1);
  mpfi_init(c_arg_tmp);
  mpfi_c_init(c_sqrt_tmp);
  mpfi_init(pp_tlnx);
  mpfi_init(pp_x_sigma);
  mpfi_init(pow_tmp);
  mpfi_init(pow_tmp1);
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
  mpfi_init(mpfi_pi_by_2);
  mpfi_div_2ui(mpfi_pi_by_2,mpfi_pi,1);
  mpfi_init(mpfi_pi_by_4);
  mpfi_div_2ui(mpfi_pi_by_4,mpfi_pi,2);
  mpfi_init(mpfi_log_2);
  mpfi_const_log2(mpfi_log_2);
  mpfi_init(mpfi_ln_pi);
  mpfi_log(mpfi_ln_pi,mpfi_pi);
}

// res=|z|^2
inline void mpfi_c_norm(mpfi_ptr res, mpfi_c_ptr z)
{
  mpfi_sqr(norm_tmp1,z->re);
  mpfi_sqr(norm_tmp,z->im);
  mpfi_add(res,norm_tmp1,norm_tmp);
}

// res=|z|
inline void mpfi_c_abs(mpfi_ptr res, mpfi_c_ptr z)
{
   mpfi_c_norm(res,z);
   mpfi_sqrt(res,res);
}


/* use the obvious method, anything cleverer loses precision */
inline void mpfi_c_mul(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_mul(m_tmp1,z2->re,z3->re);
    mpfi_mul(m_tmp2,z2->im,z3->im);
    mpfi_sub(m_tmp3,m_tmp1,m_tmp2);
    mpfi_mul(m_tmp1,z2->re,z3->im);
    mpfi_mul(m_tmp2,z2->im,z3->re);
    mpfi_add(z1->im,m_tmp1,m_tmp2);
    mpfi_swap(z1->re,m_tmp3); // for cache reasons, might actually be quicker to use mpfi_set
};

// better than mpfi_c_mul(z1,z2,z2)
inline void mpfi_c_sqr(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
  mpfi_sqr(m_tmp1,z2->re);
  mpfi_sqr(m_tmp2,z2->im);
  mpfi_mul(m_tmp3,z2->re,z2->im);
  mpfi_sub(z1->re,m_tmp1,m_tmp2);
  mpfi_mul_2ui(z1->im,m_tmp3,1);
}

/* multiply (mpfi_c_t z2) by (mpfi_t x) into (mpfi_c_t z1) safely */
inline void mpfi_c_mul_i(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_ptr x)
{
    mpfi_mul(z1->re,z2->re,x);
    mpfi_mul(z1->im,z2->im,x);
};

// scalar mult
inline void mpfi_c_mul_d(mpfi_c_ptr z1, mpfi_c_ptr z2, double x)
{
    mpfi_mul_d(z1->re,z2->re,x);
    mpfi_mul_d(z1->im,z2->im,x);
}

// and again
inline void mpfi_c_mul_ui(mpfi_c_ptr z1, mpfi_c_ptr z2, unsigned long int x)
{
    mpfi_mul_ui(z1->re,z2->re,x);
    mpfi_mul_ui(z1->im,z2->im,x);
}

/* multiply by i in situ */
inline void mpfi_c_muli(mpfi_c_ptr z)   
{
    mpfi_swap(z->re,z->im);
    mpfi_neg(z->re,z->re);
};

/* multiply (mpfi_c_t z2) by (mpfr_t x) into (mpfi_c_t z1) safely */
inline void mpfi_c_mul_fr(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfr_ptr x)
{
    mpfi_mul_fr(z1->re,z2->re,x);
    mpfi_mul_fr(z1->im,z2->im,x);
};

// scalar mult by mpz_t
inline void mpfi_c_mul_z(mpfi_c_ptr z1,mpfi_c_ptr z2, mpz_ptr x)
{
    mpfi_mul_z(z1->re,z2->re,x);
    mpfi_mul_z(z1->im,z2->im,x);
};
    

/* z1<-conjugate(z2) safely */
inline void mpfi_c_conj(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_c_set(z1,z2);
    mpfi_neg(z1->im,z1->im);
};

// z1=z2/z3
inline void mpfi_c_div(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_c_ptr z3)
{
    mpfi_sqr(d_tmp1,z3->re);
    mpfi_sqr(d_tmp2,z3->im);
    mpfi_add(d_tmp4,d_tmp1,d_tmp2);
    mpfi_c_conj(d_tmp3,z3);
    mpfi_c_mul(d_tmp5,d_tmp3,z2);
    mpfi_div(z1->re,d_tmp5->re,d_tmp4);
    mpfi_div(z1->im,d_tmp5->im,d_tmp4);
};

// scalar division by mpz
inline void mpfi_c_div_z(mpfi_c_ptr z1,mpfi_c_ptr z2, mpz_ptr x)
{
    mpfi_div_z(z1->re,z2->re,x);
    mpfi_div_z(z1->im,z2->im,x);
}

// scalar division by mpfi_t
inline void mpfi_c_div_i (mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_ptr x)
{
	mpfi_div(z1->re,z1->re,x);
	mpfi_div(z1->im,z1->im,x);

}

// scalar division by ui
inline void mpfi_c_div_ui (mpfi_c_ptr z1, mpfi_c_ptr z2, unsigned long int i)
{
    mpfi_div_ui(z1->re,z2->re,i);
    mpfi_div_ui(z1->im,z2->im,i);
};

// scalar division by double
inline void mpfi_c_div_d(mpfi_c_ptr z1, mpfi_c_ptr z2, double x)
{
  mpfi_div_d(z1->re,z2->re,x);
  mpfi_div_d(z1->im,z2->im,x);
}

// z1=exp(z2)
inline void mpfi_c_exp(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_exp(e_tmp,z2->re);
    mpfi_cos(z1->re,z2->im);
    mpfi_sin(z1->im,z2->im);
    mpfi_c_mul_i(z1,z1,e_tmp);
};

// res=x^(re_s+I im_s)
inline void mpfi_c_pow_double_to_doubles(mpfi_c_ptr res, double x, double re_s, double im_s)
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

// res=x^y in reals
inline void mpfi_pow(mpfi_ptr res, mpfi_ptr x, mpfi_ptr y)
{
  mpfi_log(pow_tmp,x);
  mpfi_mul(pow_tmp1,pow_tmp,y);
  mpfi_exp(res,pow_tmp1);
}

// res=x^z Complex=Real^Complex
inline void mpfi_c_pow_i_c (mpfi_c_ptr res, mpfi_ptr x, mpfi_c_ptr z)
{
  mpfi_pow(pp_x_sigma,x,z->re);
  mpfi_log(pp_tlnx,x);
  mpfi_mul(pp_tlnx,pp_tlnx,z->im);
  mpfi_cos(res->re,pp_tlnx);
  mpfi_sin(res->im,pp_tlnx);
  mpfi_c_mul_i(res,res,pp_x_sigma);
}

// this appear to have got defined in mpfi somewhere along the way
/*
inline void mpfi_atan2 (mpfi_ptr res, mpfi_ptr y, mpfi_ptr x)
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

 // res=log(z)
inline void mpfi_c_log (mpfi_c_ptr res, mpfi_c_ptr z)
{
  mpfi_c_norm(c_log_tmp,z);
  mpfi_log(c_log_tmp1,c_log_tmp);
  mpfi_div_ui(c_log_tmp,c_log_tmp1,2);
  mpfi_atan2(res->im,z->im,z->re);
  mpfi_set(res->re,c_log_tmp);
}

// aka mpfi_c_abs
#define mpfi_c_mod(a,b) mpfi_c_abs(a,b)

// arg(z)
inline void mpfi_c_arg(mpfi_ptr res, mpfi_c_ptr z)
{
  mpfi_div(c_arg_tmp,z->im,z->re);
  mpfi_atan(res,c_arg_tmp);
}

inline int mpfi_contains_zero(mpfi_ptr x)
{
  return(mpfi_cmp_d(x,0.0)==0);
}

inline int mpfi_c_contains_zero(mpfi_c_ptr z)
{
  return(mpfi_contains_zero(z->re)&&mpfi_contains_zero(z->im));
}

//z=0+0i
inline void mpfi_c_zero(mpfi_c_ptr z)
{
   mpfi_set_ui(z->re,0);
   mpfi_set_ui(z->im,0);
}

//z=1+0i
inline void mpfi_c_one(mpfi_c_ptr z)
{
  mpfi_set_ui(z->re,1);
  mpfi_set_ui(z->im,0);
}

// returns root in right half plane
inline void mpfi_c_sqrt(mpfi_c_ptr res, mpfi_c_ptr z)
{
   mpfi_c_log(c_sqrt_tmp,z);
   mpfi_div_2ui(c_sqrt_tmp->re,c_sqrt_tmp->re,1);
   mpfi_div_2ui(c_sqrt_tmp->im,c_sqrt_tmp->im,1);
   mpfi_c_exp(res,c_sqrt_tmp);
   if(mpfi_is_neg(res->re))
     mpfi_c_neg(res,res);
}

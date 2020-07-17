

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

void mpfi_print(mpfi_ptr x)
{
    mpfi_out_str(stdout,10,0,x);
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

void mpfi_c_set_ui(mpfi_c_ptr z, unsigned int re, unsigned int im)
{
  mpfi_set_ui(z->re,re);
  mpfi_set_ui(z->im,im);
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

void mpfi_c_sub_re(mpfi_c_ptr z1, mpfi_c_ptr z2, mpfi_ptr x)
{
    mpfi_sub(z1->re,z2->re,x);
};

void mpfi_c_dec_ui(mpfi_c_ptr z, unsigned int i)
{
  mpfi_sub_ui(z->re,z->re,i);
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

mpfi_c_t d_tmp3,F_tmp1,F_tmp2;
mpfi_t m_tmp1,m_tmp2,m_tmp3,a_tmp,d_tmp1,d_tmp2,e_tmp;   

void mpfi_c_setup(int prec)
{
/* now initialise all global variables.
   we don't bother clearing them down afterwards, just exit. */
    mpfr_set_default_prec(prec);
    mpfi_init(m_tmp1);
    mpfi_init(m_tmp2);
    mpfi_init(m_tmp3);
    mpfi_init(d_tmp1);
    mpfi_init(d_tmp2);
    mpfi_init(a_tmp);
    mpfi_c_init(d_tmp3);
    mpfi_init(e_tmp);
    mpfi_c_init(F_tmp1);
    mpfi_c_init(F_tmp2);

};

void mpfi_c_abs(mpfi_ptr res, mpfi_c_ptr z)
{
   mpfi_sqr(res,z->re);
   mpfi_sqr(a_tmp,z->im);
   mpfi_add(res,res,a_tmp);
   mpfi_sqrt(res,res);
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

void mpfi_c_div_ui (mpfi_c_ptr z1, mpfi_c_ptr z2, unsigned long int i)
{
    mpfi_div_ui(z1->re,z2->re,i);
    mpfi_div_ui(z1->im,z2->im,i);
};

void mpfi_c_exp(mpfi_c_ptr z1, mpfi_c_ptr z2)
{
    mpfi_exp(e_tmp,z2->re);
//    mpfi_sin_cos(z1->im,z1->re,z2->im);   would be nice, but does not exist
    mpfi_cos(z1->re,z2->im);
    mpfi_mul(z1->re,z1->re,e_tmp);
    mpfi_sin(z1->im,z2->im);
    mpfi_mul(z1->im,z1->im,e_tmp);
};


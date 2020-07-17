#ifndef MPFI_C
#define MPFI_C

// First the mpfi_c stuff

/* define a complex version of an MPFI interval */
typedef struct{
    mpfi_t re;
    mpfi_t im;
} _mpfi_c_struct;

// structs to make parameter passing easier
typedef _mpfi_c_struct mpfi_c_t[1];
typedef _mpfi_c_struct *mpfi_c_ptr;

mpfi_t mpfi_pi,mpfi_2_pi;
/* initialisation */
inline void mpfi_c_init(mpfi_c_ptr);

/* clearing */
inline void mpfi_c_clear(mpfi_c_ptr);

/* efficient swapping */
inline void mpfi_c_swap(mpfi_c_ptr,mpfi_c_ptr);

/* simple printing */
inline void mpfi_c_print(mpfi_c_ptr);

inline void mpfi_c_print_str(const char *, mpfi_c_ptr);

inline void mpfi_c_printn(mpfi_c_ptr,int);

inline void mpfi_print(mpfi_ptr);

inline void mpfi_print_str(const char *, mpfi_ptr);

inline void mpfi_printn(mpfi_ptr, int);

// setting from various data types
inline void mpfi_c_set(mpfi_c_ptr, mpfi_c_ptr);

inline void mpfi_c_set_d(mpfi_c_ptr, double, double);

inline void mpfi_c_set_ui(mpfi_c_ptr, unsigned long int, unsigned long int);

// just set the real part
inline void mpfi_c_set_re(mpfi_c_ptr,mpfi_ptr);

inline void mpfi_c_set_re_d(mpfi_c_ptr,double);

// just set the imaginary part
inline void mpfi_c_set_im(mpfi_c_ptr, mpfi_ptr);

// z1=z2+z3
inline void mpfi_c_add(mpfi_c_ptr, mpfi_c_ptr, mpfi_c_ptr);

// z1=z2-z3
inline void mpfi_c_sub(mpfi_c_ptr, mpfi_c_ptr, mpfi_c_ptr);

// z1=-z2
inline void mpfi_c_neg(mpfi_c_ptr,mpfi_c_ptr);

inline void mpfi_c_setup(unsigned long int);

// res=|z|^2
inline void mpfi_c_norm(mpfi_ptr, mpfi_c_ptr);

// res=|z|
inline void mpfi_c_abs(mpfi_ptr, mpfi_c_ptr);


/* use the obvious method, anything cleverer loses precision */
inline void mpfi_c_mul(mpfi_c_ptr, mpfi_c_ptr, mpfi_c_ptr);

// better than mpfi_c_mul(z1,z2,z2)
inline void mpfi_c_sqr(mpfi_c_ptr, mpfi_c_ptr);

/* multiply (mpfi_c_t z2) by (mpfi_t x) into (mpfi_c_t z1) safely */
inline void mpfi_c_mul_i(mpfi_c_ptr, mpfi_c_ptr, mpfi_ptr);

// scalar mult
inline void mpfi_c_mul_d(mpfi_c_ptr, mpfi_c_ptr, double);

// and again
inline void mpfi_c_mul_ui(mpfi_c_ptr, mpfi_c_ptr, unsigned long int);

/* multiply by i in situ */
inline void mpfi_c_muli(mpfi_c_ptr);

/* multiply (mpfi_c_t z2) by (mpfr_t x) into (mpfi_c_t z1) safely */
inline void mpfi_c_mul_fr(mpfi_c_ptr, mpfi_c_ptr, mpfr_ptr);

// scalar mult by mpz_t
inline void mpfi_c_mul_z(mpfi_c_ptr,mpfi_c_ptr, mpz_ptr);
    

/* z1<-conjugate(z2) */
inline void mpfi_c_conj(mpfi_c_ptr, mpfi_c_ptr);

// z1=z2/z3
inline void mpfi_c_div(mpfi_c_ptr, mpfi_c_ptr, mpfi_c_ptr);

// scalar division by mpz
inline void mpfi_c_div_z(mpfi_c_ptr,mpfi_c_ptr, mpz_ptr);

// scalar division by mpfi_t
inline void mpfi_c_div_i (mpfi_c_ptr, mpfi_c_ptr, mpfi_ptr);

// scalar division by ui
inline void mpfi_c_div_ui (mpfi_c_ptr, mpfi_c_ptr, unsigned long int);

// scalar division by double
inline void mpfi_c_div_d(mpfi_c_ptr, mpfi_c_ptr, double);

// z1=exp(z2)
inline void mpfi_c_exp(mpfi_c_ptr, mpfi_c_ptr);

// res=x^(re_s+I im_s)
inline void mpfi_c_pow_double_to_doubles(mpfi_c_ptr, double, double, double);

// res=x^y in reals
inline void mpfi_pow(mpfi_ptr, mpfi_ptr, mpfi_ptr);

// res=x^z Complex=Real^Complex
inline void mpfi_c_pow_i_c (mpfi_c_ptr, mpfi_ptr, mpfi_c_ptr);

 // res=log(z)
inline void mpfi_c_log (mpfi_c_ptr, mpfi_c_ptr);

// aka mpfi_c_abs
#define mpfi_c_mod(a,b) mpfi_c_abs(a,b)

// arg(z)
inline void mpfi_c_arg(mpfi_ptr, mpfi_c_ptr);

inline int mpfi_contains_zero(mpfi_ptr);

inline int mpfi_c_contains_zero(mpfi_c_ptr);

//z=0+0i
inline void mpfi_c_zero(mpfi_c_ptr);

//z=1+0i
inline void mpfi_c_one(mpfi_c_ptr);

// returns root in right half plane
inline void mpfi_c_sqrt(mpfi_c_ptr, mpfi_c_ptr);
#endif

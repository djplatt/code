#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
//#include "mpfi.h"
//#include "mpfi_io.h"
//#include "../includes/mpfi_c.h"
//#include "../includes/mpfi_fft.h"
//#include "../includes/fft_defs.h"

typedef struct {mpfr_t left;mpfr_t right;} mpfi_t;

typedef _mpfi_struct mpfi_t[1];
typedef _mpfi_struct *mpfi_ptr;

mpfr_t mpfr_tmp;
mpfi_t mpfi_tmp;

#define mpfi_init(x) {mpfr_init(x->left);mpfr_init(x->right);}

#define mpfi_clear(x) {mpfr_clear(x->left);mpfr_clear(x->right);}

#define mpfi_set(x,y) {mpfr_set(x->left,y->left,GMP_RNDD);mpfr_set(x->right,y->right,GMP_RNDU);}

#define mpfi_set_ui(x,y) {mpfr_set_ui(x->left,y,GMP_RNDD);mpfr_set_ui(x->right,y,GMP_RNDU);}

#define mpfi_add(x,y,z) {mpfr_add(x->left,y->left,z->left,GMP_RNDD);mpfr_add(x->right,y->right,z->right,GMP_RNDU);}

#define mpfi_neg(x,y) {mpfr_neg(mpfr_tmp,y->right,GMP_RNDD);mpfr_neg(x->right,y->left,GMP_RNDU);mpfr_set(x->left,mpfr_tmp,GMP_RNDD);}

#define mpfi_sub(x,y,z) {mpfr_sub(mpfi_tmp->left,y->left,z->right,GMP_RNDD);mpfr_sub(mpfi_tmp->right,y->right,z->left,GMP_RNDU);mpfi_set(x,mpfi_tmp);}

inline void mpfi_mul(mpfi_ptr x, mpfi_ptr y, mpfi_ptr z) 
{
  int syl=mpfr_sgn(y->left),szl=mpfr_sgn(z->left);
  if(syl>=0) // y > 0
    {
      if(szl>=0) // (y >= 0) && (z >= 0)
	{
	  mpfr_mul(x->left,y->left,z->left,GMP_RNDD);
	  mpfr_mul(x->right,y->right,z->right,GMP_RNDU);
	  return;
	}
      else // (y >= 0) &&((z <=0) or (0 in z))
	{
	  mpfr_mul(mpfr_tmp,y->left,z->right,GMP_RNDD);
	  mpfr_mul(x->right,y->right,z->left,GMP_RNDU);
	  mpfr_set(x->left,mpfr_tmp,GMP_RNDD);
	  return;
	}
    }
  if(szl>=0) // (z>=0) && ((y <=0 ) or (0 in y))
    {	      
      mpfr_mul(mpfr_tmp,y->left,z->right,GMP_RNDD);
      mpfr_mul(x->right,y->right,z->left,GMP_RNDU);
      mpfr_set(x->left,mpfr_tmp,GMP_RNDD);
      return;
    }




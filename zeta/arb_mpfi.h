#include "arb.h"
#include "acb.h"
int64_t PREC;
#define mpfi_t arb_t
#define mpfi_ptr arb_ptr
#define mpfi_init(x) arb_init(x)
#define mpfi_clear(x) arb_clear(x)

#define mpfi_is_neg(x) arb_is_nonpositive(x)
#define mpfi_is_pos(x) arb_is_nonnegative(x)

#define mpfi_set_ui(x,y) arb_set_ui(x,y)
#define mpfi_set_d(x,y) arb_set_d(x,y)
#define mpfi_neg(x,y) arb_neg(x,y)
#define mpfi_add_ui(x,y,z) arb_add_ui(x,y,z,PREC)

#define mpfi_mul_2ui(x,y,e) arb_mul_2exp_si(x,y,e)
#define mpfi_div_2ui(x,y,e) arb_mul_2exp_si(x,y,-((int64_t) e))

#define mpfi_c_t acb_t
#define mpfi_c_ptr acb_ptr
#define mpfi_c_init(x) acb_init(x)
#define mpfi_c_clear(x) acb_clear(x)

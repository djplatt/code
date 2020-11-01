#ifndef WIN_ZETA
#define WIN_ZETA

#include "inttypes.h"
#include "acb.h"


// defined so that POS&NEG == 0 but no other combination
typedef int sign_t;
#define POS (1)
#define NEG (2)
#define UNK (3)

typedef int dir_t;
#define UP (1)
#define DOWN (2)

sign_t arb_sign(arb_t);

dir_t arb_dir(arb_ptr, arb_ptr);

sign_t acb_sign(acb_t);

dir_t acb_dir(acb_ptr, acb_ptr);

double arb_get_d(arb_t);

void arb_mul_d(arb_t, arb_t, double, int64_t);

void arb_div_d(arb_t, arb_t, double, int64_t);

void arb_add_d(arb_t, arb_t, double, int64_t);

void arb_sub_d(arb_t, arb_t, double, int64_t);

int arb_cmp(arb_t, arb_t);

#endif 

#ifndef INTER
#define INTER
#include "acb.h"
#include "parameters.h"

// compute exp(-t^2/2H^2)
void inter_gaussian(arb_ptr, arb_ptr, int64_t);

// compute sinc(2*B*Pi*t) into ress, cos(2*B*Pi*t) into resc
// first time we call this (with sincp true) we compute sin and cos
// from then on (sincp false) just swap signs each time
void inter_sinc_cos(arb_ptr, arb_ptr, arb_ptr, int64_t);

// t is distance from t0 to t implied by f_ptr
// if fd_res is NULL, don't calc differential
void inter_point(arb_ptr, arb_ptr, acb_t *, arb_ptr, int, int64_t);

// t0 points into f_vec
// let t=t0+(t_ptr-N/2)*one_over_A
// returns f(t) in f_res, f'(t) in fd_res
// return f_vec[t0]->re if t_ptr is integral 
void arb_inter_t(arb_ptr, arb_ptr, acb_t *, arb_ptr, int64_t);

void acb_inter_t(arb_ptr, acb_t *, arb_ptr, int64_t);
#endif

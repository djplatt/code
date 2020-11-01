#ifndef TURING
#define TURING
#include "inttypes.h"
#include "stdbool.h"

// we have found a stationary pt between left, left+1 and left+2
// we don't assume RN, so go find the zeros
bool resolve_stat_point(arb_ptr, arb_ptr, arb_ptr, arb_ptr, arb_ptr, arb_ptr, int, sign_t, acb_t, int64_t);

// count zeros (using stationary points) between start and end
int zeros_st(acb_t, int, int, int64_t);

// int Im loggamma(1/4+it/2)
void im_int1(arb_t, arb_t, int64_t);

void im_int(arb_t, arb_t, arb_t, int64_t);

// true if there is a stationary point l-m-r
// we have already checked for sign changes
// On RH we either have two zeros between l and m or between m and r
bool stat_pt(acb_ptr, acb_ptr, acb_ptr);

// assume that if a zero lies between t0 and t1 then it is at t0
double Nleft_int(long int, long int, acb_t, double, int64_t);

// assume that if a zero lies between t0 and t1 then it is at t1
double Nright_int(long int, long int, acb_t, double, int64_t);

// returns int_t0^{t0+h} S(t) dt
// t0=t0_ptr*delta
// t0>168*Pi
// Trudgian Thesis Theorem 5.2.2
// < 0.059 log(t) + 2.067
void St_int(arb_t, arb_t, int64_t);

// the log term integrated
void ln_term(arb_t, arb_t, arb_t, arb_t, int64_t);

// the maximum number of zeros <=a based on Turing region [a,b]
long int turing_max(acb_t, long int, double, long int, double, int64_t);

long int turing_min(acb_t, long int, double, long int, double, int64_t);

// Use Turing's method to estimate and then check number of zeros in [a,b]
// last_max is zeros to a (if known, 0 otherwise)
// on success returns the number of zeros of zeta up to b
// on failure return 0
// a_ptr,b_ptr are pointers to a,b respectively in f_vec
long int turing(acb_t *, long int, double, long int, double, long int, int64_t);

#endif

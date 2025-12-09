#include "flint/acb.h"
#include "inttypes.h"
void arb_maxd2(arb_ptr res, 
	       void (*f1)(acb_t,const acb_t,int64_t), 
	       void (*f2)(acb_t,const acb_t,int64_t), 
	       const arb_ptr low, const arb_ptr hi,
	       int64_t prec,
	       uint64_t steps);
void arb_maxd(arb_ptr res, 
	       void (*f1)(acb_t,const acb_t,int64_t), 
	       const arb_ptr low, const arb_ptr hi,
	       int64_t prec,
	       uint64_t steps);

void comp_error(arb_t err, const arb_t maxd, int64_t n, int64_t prec);

void molin_int(arb_ptr res, int64_t n, 
	       void (*f1)(arb_t,const arb_t,int64_t), 
	       const arb_ptr maxd,
	       const arb_ptr low, const arb_t hi,
	       int64_t prec);

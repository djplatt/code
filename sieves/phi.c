#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../code/includes/mpfi_c.h"
#include "math.h"
#include "inttypes.h"
#include "../code/includes/pi_x.h"

//#define DEBUG
//#define TIME

#ifdef DEBUG
#define ifdebug(x) x
#else
#define ifdebug(x)
#endif

#ifdef TIME
#include "time.h"
#endif

#define TRUE (0==0)
#define FALSE (1==0)
#define FAILURE (1)
#define bool int

// convert a 128 bit unsigned int into mpfi_t
void mpfi_set_128(mpfi_ptr res, __uint128_t x)
{
  __uint64_t tmp=(__uint64_t) x; // just gets bottom 64 bits
  ifdebug(printf("bottom 64 bits = %lu\n",tmp));
  mpfi_set_ui(res,x>>64); // gets top 64 bits
  mpfi_mul_2ui(res,res,64);
  mpfi_add_ui(res,res,tmp);
}

mpfi_t lambdas[7],tlogs[3],t0s[4],root2,root2pi,texp,ta3,ttmp,ttmp1;
mpfi_t ptmp;

void init_taylors(double lambda)
{
  ptype i;
  mpfi_c_setup(PREC);
  for(i=0;i<7;i++)
    mpfi_init(lambdas[i]);
  mpfi_set_d(lambdas[0],lambda);
  for(i=1;i<7;i++)
    mpfi_mul(lambdas[i],lambdas[i-1],lambdas[0]);
  for(i=0;i<3;i++)
    {
      mpfi_init(tlogs[i]);
      mpfi_init(t0s[i]);
    }
  mpfi_init(t0s[3]);
  mpfi_init(root2);
  mpfi_init(root2pi);
  mpfi_set_ui(root2,2);
  mpfi_sqrt(root2,root2);
  mpfi_mul(root2pi,mpfi_sqrt_pi,root2);
  mpfi_init(texp);
  mpfi_init(ta3);
  mpfi_init(ttmp);
  mpfi_init(ttmp1);
  mpfi_init(ptmp);
}


	      
// use when t0<x
// phi(t)=1-1/2 erfc(log(t/x)/(sqrt(2)*lambda))
//       =1/2-erf(log(t/x)/(sqrt(2)*lambda))
void taylor_coeff_left(mpfi_ptr a0, mpfi_ptr a1, mpfi_ptr a2, 
		  mpfi_ptr err, mpfi_ptr x, mpfi_ptr t0, ptype xi)
{
  ptype i;
  mpfi_set(t0s[0],t0); // t0
  mpfi_mul(t0s[1],t0s[0],t0s[0]); // t0^2
  mpfi_mul(t0s[2],t0s[1],t0s[0]); // t0^3
  mpfi_mul(t0s[3],t0s[1],t0s[1]); // t0^4

  mpfi_div(tlogs[0],t0s[0],x);
  mpfi_log(tlogs[0],tlogs[0]); // t0<x => log(t0/x)<0
  mpfi_sqr(tlogs[1],tlogs[0]); // log^2
  mpfi_mul(tlogs[2],tlogs[1],tlogs[0]); // log^3
  mpfi_sqr(tlogs[3],tlogs[1]); // log^4

  mpfi_div(a0,tlogs[0],root2);
  mpfi_div(a0,a0,lambdas[0]);
  mpfi_erf(a0,a0);
  mpfi_div_2ui(a0,a0,1);
  mpfi_add_d(a0,a0,0.5); // 1/2-erf(log(t0/x)/sqrt(2)/lambda)

  mpfi_div(texp,tlogs[1],lambdas[1]); // log^2/lam^2
  mpfi_mul_d(texp,texp,-0.5);
  mpfi_exp(texp,texp); // exp(-log^2/2/lam^2) 
#ifdef DEBUG
  mpfi_print_str("exp()=",texp);
#endif
  mpfi_mul(ttmp,root2pi,t0s[0]); // sqrt(2 pi)*t0
  mpfi_mul(ttmp,ttmp,lambdas[0]); // sqrt(2 pi)*t0*lam

  ifdebug(mpfi_print_str("denom=",ttmp));
  mpfi_div(a1,texp,ttmp); // -exp(-log^2/2/lam^2)/(2 sqrt(2 pi)*t0*lam)

  mpfi_add(a2,lambdas[1],tlogs[0]);
  mpfi_mul(a2,a2,texp); // exp()*(log+lam^2)

  mpfi_mul(ttmp,ttmp,lambdas[1]);
  mpfi_mul(ttmp,ttmp,t0s[0]); // sqrt(2 pi)*t0^2*lam^3
  mpfi_mul_ui(ttmp,ttmp,2); // 2*sqrt....
  mpfi_div(a2,a2,ttmp);
  mpfi_neg(a2,a2);

  mpfi_mul_ui(ta3,tlogs[0],3);
  mpfi_sub_ui(ta3,ta3,1);
  mpfi_mul(ta3,ta3,lambdas[1]); // lam^2(3*log-1)
  ifdebug(mpfi_print_str("(3*log-1)lam^2=",ta3));
  mpfi_add(ta3,ta3,lambdas[3]);
  mpfi_add(ta3,ta3,lambdas[3]);
  ifdebug(mpfi_print_str("+2lam^4=",ta3));
  mpfi_add(ta3,ta3,tlogs[1]);
  ifdebug(mpfi_print_str("(2*lam^4...)=",ta3));

  mpfi_mul(ta3,ta3,texp);
  ifdebug(mpfi_print_str("a3 num=",ta3));

  mpfi_mul(ttmp,ttmp,lambdas[1]);
  mpfi_mul(ttmp,ttmp,t0s[0]);
  mpfi_mul_ui(ttmp,ttmp,3);
  mpfi_div(ta3,ta3,ttmp);
  ifdebug(mpfi_print_str("a3=",ta3));
  mpfi_mul_d(err,ta3,0.25);
  mpfi_mul_ui(err,err,xi);
  mpfi_mul_ui(err,err,xi);
  mpfi_mul_ui(ttmp1,err,3);
  mpfi_add(a1,a1,ttmp1);
  mpfi_mul_ui(err,ttmp1,xi);
  // err is now the residual error from the a3 term
  // need to add the error for the a4 and subs. terms

  mpfi_sub_ui(tlogs[0],t0s[0],xi);
  mpfi_div(tlogs[0],tlogs[0],x);
  mpfi_add_ui(ttmp,t0s[0],xi);
  mpfi_div(ttmp,ttmp,x);
  mpfi_put(tlogs[0],ttmp);
  mpfi_log(tlogs[0],tlogs[0]);
  ifdebug(mpfi_print_str("log([t0-xi,t0+xi]/x)=",tlogs[0]));
  mpfi_sqr(tlogs[1],tlogs[0]);
  mpfi_mul(tlogs[2],tlogs[1],tlogs[0]);
  mpfi_div(texp,tlogs[1],lambdas[1]);
  mpfi_mul_d(texp,texp,-0.5);
  mpfi_exp(texp,texp);
  ifdebug(mpfi_print_str("exp term =",texp));
  mpfi_set_ui(ttmp,xi);
  mpfi_sqr(ttmp,ttmp);
  mpfi_sqr(ttmp,ttmp);
  mpfi_mul(texp,texp,ttmp); // exp(-..)*xi^4
  mpfi_sub_ui(ttmp,t0s[0],xi);
  mpfi_sqr(ttmp,ttmp);
  mpfi_sqr(ttmp,ttmp);
  mpfi_mul(ttmp,ttmp,root2pi);
  mpfi_mul_ui(ttmp,ttmp,24);
  mpfi_mul(ttmp,ttmp,lambdas[6]); // 24sqrt(2pi)*lam^7*(t0-xi)^4
  mpfi_div(texp,texp,ttmp);
  mpfi_mul_ui(ta3,tlogs[0],11);
  mpfi_sub_ui(ttmp,ta3,6);
  mpfi_mul(ta3,ttmp,lambdas[3]); // (11log()-6)lam^4
  mpfi_mul_ui(ttmp,tlogs[0],3);
  mpfi_sub(ttmp,tlogs[1],ttmp);
  mpfi_mul_ui(ttmp,ttmp,6);
  mpfi_mul(ttmp,ttmp,lambdas[1]); // 6(log^2()-3log)lam^2
  mpfi_add(ta3,ta3,ttmp);
  mpfi_add(ta3,ta3,tlogs[2]); // + log^3()
  mpfi_mul_ui(ttmp,lambdas[5],6); // 6lam^6
  mpfi_add(ta3,ta3,ttmp);
  mpfi_mul(ttmp,ta3,texp);
  ifdebug(mpfi_print_str("Additional Taylor error=",ttmp));
  mpfi_add(err,err,ttmp);
  mpfi_neg(ttmp,err);
  mpfi_put(err,ttmp);
  ifdebug(mpfi_print_str("Total Taylor Error=",err));
}

void phi(mpfi_ptr res, ptype s1, __int64_t st, bigint st2, mpfi_ptr a0, mpfi_ptr a1, mpfi_ptr a2, mpfi_ptr err)
{
  __uint64_t tmp;
  if(st<0)
    {
      tmp=-st;
      mpfi_mul_ui(res,a1,tmp); // a1 * sigma (t-t0)
      mpfi_neg(res,res);
    }
  else
    mpfi_mul_ui(res,a1,st); // a1 * sigma (t-t0)

  ifdebug({printf("s1=%lu\nst=%ld\nst2=",s1,st);print_bigint(st2);printf("\n");});
  ifdebug(mpfi_print_str("a1*sigma (t-t0)",res));
  mpfi_set_128(ptmp,st2); // sigma (t-t0)^2
  ifdebug(mpfi_print_str("sigma (t-t0)^2=",ptmp));
  mpfi_mul(ptmp,a2,ptmp); // a2* ...
  mpfi_add(res,res,ptmp);
  mpfi_mul_ui(ptmp,a0,s1);   // sigma 1
  mpfi_add(res,res,ptmp);
  mpfi_mul_ui(ptmp,err,s1);  // err * sigma 1
  mpfi_add(res,res,ptmp);
}
  
int main()
{

  ptype i;
  mpfi_t x,t0,a0,a1,a2,err,res;
  bignum xb=calc_x();

  init_taylors(LAMBDA);
  mpfi_init(x);
  mpfi_init(t0);
  mpfi_init(a0);mpfi_init(a1);mpfi_init(a2);mpfi_init(err);
  mpfi_init(res);
  mpfi_set_128(x,xb);

  for(i=0;i<1000000;i++)
    {
      taylor_coeff_left(a0,a1,a2,err,x,t0,XI);
      //mpfi_print_str("a0=",a0);
      //mpfi_print_str("a1=",a1);
      //mpfi_print_str("a2=",a2);
      //mpfi_print_str("err=",err);

      phi(res,1,-XI,(bigint) XI * (bigint) XI, a0,a1,a2,err);
      //mpfi_print_str("phi(t0-xi)=",res);
      //phi(res,1,XI,(bigint) XI * (bigint) XI, a0,a1,a2,err);
      //mpfi_print_str("phi(t0+xi)=",res);
    }
  return(0);
}

  


  

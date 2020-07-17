#include "arb.h"
#include <primesieve.hpp>
#include "inttypes.h"


// ms = -s
void g(arb_t res, uint64_t p, arb_t ms, int64_t prec)
{
  static bool init=false;
  static arb_t ap,ps,tmp1,tmp2,tmp3,tmp4;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(tmp4);
      arb_init(ps);
      arb_init(ap);
    }
  arb_set_ui(ap,p);
  arb_pow(ps,ap,ms,prec);
  arb_set_ui(tmp1,4);
  arb_div_ui(tmp4,tmp1,p,prec);
  arb_div_ui(tmp2,tmp4,(p-2),prec);
  arb_mul(tmp1,tmp2,ps,prec); // 4/p*(p-2) p^-s
  //arb_printd(tmp1,10);printf("\n");
  arb_set_ui(tmp2,3*p+2);
  arb_div_ui(tmp3,tmp2,p,prec);
  arb_div_ui(tmp4,tmp3,p,prec);
  arb_div_ui(tmp3,tmp4,p-2,prec);
  arb_mul(tmp2,tmp3,ps,prec);
  arb_mul(tmp3,tmp2,ps,prec); // (3p+2)/p^2(p-2) p^-2s
  //arb_printd(tmp3,10);printf("\n");
  arb_add(tmp2,tmp1,tmp3,prec);
  arb_set_ui(tmp1,2);
  arb_div_ui(tmp3,tmp1,p,prec);
  arb_div_ui(tmp4,tmp3,p,prec);
  arb_div_ui(tmp3,tmp4,p-2,prec);
  arb_mul(tmp1,tmp3,ps,prec);
  arb_mul(tmp3,tmp1,ps,prec);
  arb_mul(tmp1,tmp3,ps,prec);
  //arb_printd(tmp1,10);printf("\n");
  arb_add(tmp3,tmp1,tmp2,prec);
  //arb_printd(tmp3,10);printf("\n");
  arb_log1p(res,tmp3,prec);
}

arb_t tmp,sum,ms;
int64_t prec;

void callback(uint64_t p)
{
  g(tmp,p,ms,prec);
  //printf("g(%lu)=",p);arb_printd(tmp,10);printf("\n");
  arb_add(sum,sum,tmp,prec);
}

int main(int argc, char **argv)
{
  if(argc!=5)
    {
      printf("Usage:- %s <P> <s num> <s denom> <prec>\n",argv[0]);
      return 0;
    }
  uint64_t P=atol(argv[1]);
  prec=atol(argv[4]);
  arb_init(sum);
  arb_init(tmp);
  arb_t tmp1,tmp2,tmp3;
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);
  arb_init(ms);
  arb_set_ui(ms,atol(argv[2]));
  arb_div_ui(ms,ms,atol(argv[3]),prec);
  printf("s=");arb_printd(ms,20);printf("\n");
  arb_set_ui(tmp,2);
  arb_pow(tmp1,tmp,ms,prec); // 2^-s
  arb_mul(tmp,tmp1,tmp1,prec); // 2^-2s
  arb_mul_ui(tmp2,tmp,3,prec); // 3.2^-2s
  arb_mul_2exp_si(tmp2,tmp2,-2); // 3/4 2^-2s
  arb_mul(tmp3,tmp,tmp1,prec); // 2^-3s
  arb_mul_2exp_si(tmp3,tmp3,-2); // 1/4 2^-3s
  arb_add(tmp,tmp3,tmp2,prec);
  arb_log1p(sum,tmp,prec);

  primesieve::callback_primes(3,P,callback);

  arb_neg(ms,ms);
  printf("log H_%ld(",P);arb_printd(ms,10);printf(") = ");arb_printd(sum,20);printf("\n");

  arb_exp(tmp,sum,prec);
  printf("H_%ld(",P);arb_printd(ms,10);printf(") = ");arb_printd(tmp,20);printf("\n");

  return 0;
}

// sum from stdin a b e
// where x =[a,b]*2^e
// a,b,e fmpz_t (signed decimal integers)
#include "stdio.h"
#include "stdlib.h"
#include "acb.h"
#include "inttypes.h"

uint64_t arb_sum(arb_ptr res, int64_t prec, fmpz_t a, fmpz_t b, fmpz_t e)
{
  uint64_t count=0;
  arf_t l,r;
  arf_init(l);
  arf_init(r);
  arb_t tmp;
  arb_init(tmp);
  arb_zero(tmp);
  uint64_t st;
  while(fscanf(stdin,"%lu ",&st)==1)
    {
      if((fmpz_read(a)<=0)||(fmpz_read(b)<=0)||(fmpz_read(e)<=0))
	{
	  arb_clear(tmp);
	  arf_clear(l);
	  arf_clear(r);
	  printf("Error reading a,b,e to go with a in artur_arb_sum. Exiting.\n");
	  exit(0);
	}
      count++;
      arf_set_fmpz(l,a);
      arf_set_fmpz(r,b);
      arf_mul_2exp_fmpz(l,l,e);
      arf_mul_2exp_fmpz(r,r,e);
      arb_set_interval_arf(tmp,l,r,prec);
      arb_add(res,res,tmp,prec);
      printf("%lu ",st);arb_printd(res,30);printf("\n");
    }
  arb_get_interval_fmpz_2exp(a,b,e,res);

  arb_clear(tmp);
  arf_clear(l);
  arf_clear(r);
  return(count);
}

int main()
{
  fmpz_t a,b,e;
  arb_t res;
  fmpz_init(a);
  fmpz_init(b);
  fmpz_init(e);
  arb_init(res);
  arb_sum(res,200,a,b,e);
  return 0;
}

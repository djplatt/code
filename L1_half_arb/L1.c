#include "inttypes.h"
#include "acb_dirichlet.h"
#define MAX_Q (10)

uint64_t gcd (uint64_t a, uint64_t b)
/* Euclid algorithm gcd */
{
  uint64_t c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
};

int co_prime(uint64_t a, uint64_t b)
{
  return(gcd(a,b)==1);
};


int main()
{
  int64_t prec=200;
  dirichlet_group_t G;
  acb_dirichlet_hurwitz_precomp_t pre;
  acb_t tmp;
  acb_init(tmp);
  arb_set_d(acb_realref(tmp),1.0);

  uint64_t A,K,N;
  acb_dirichlet_hurwitz_precomp_choose_param(&A,&K,&N,tmp,10000.0,prec);
  acb_dirichlet_hurwitz_precomp_init(pre,tmp,1,A,K,N,prec);
  acb_clear(tmp);

  
  acb_t w[MAX_Q],v[MAX_Q];
  
  for(uint64_t i=0;i<MAX_Q;i++)
    {
      acb_init(w[i]);
      acb_init(v[i]);
    }

  /*
  acb_set_d(w[0],0.5);
  acb_set_d(w[1],0.25);
  acb_dirichlet_hurwitz(w[2],w[0],w[1],prec);
  printf("With hurwitz directly ");acb_printd(w[2],20);

  acb_dirichlet_hurwitz_precomp_eval(w[2],pre,1,4,prec);
  printf("\nWith pre-hurwitz ");acb_printd(w[2],20);
  return 0;
  */

  dirichlet_char_t ch;

  
  for(uint64_t q=3;q<=MAX_Q;q++)
    {
      if((q%4)==2) continue;
      dirichlet_group_init(G,q);


  
      for(uint64_t i=1;i<q;i++)
	if(co_prime(i,q))
	  acb_dirichlet_hurwitz_precomp_eval(v[i],pre,i,q,prec);

      acb_dirichlet_dft((acb_ptr) w,(acb_srcptr) v,G,prec);

      dirichlet_char_init(ch,G);
      dirichlet_char_first_primitive(ch,G);
      while(1==1)
	{
	  uint64_t i=dirichlet_char_exp(G,ch);
	  acb_div_ui(w[i],w[i],q,prec);
	  printf("%lu %lu: ",q,i);
	  acb_printd(w[i],20);
	  printf("\n");
	  if(dirichlet_char_next_primitive(ch,G)<0)
	    break;
	}
      dirichlet_char_clear(ch); 
      dirichlet_group_clear(G);
    }

  for(uint64_t i=0;i<MAX_Q;i++)
    {
      acb_clear(w[i]);
      acb_clear(v[i]);
    }

  acb_dirichlet_hurwitz_precomp_clear(pre);  
  
  return 0;
}

#include "inttypes.h"
#include <stdlib.h>
#include "flint/acb_dirichlet.h"

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


int main(int argc, char **argv)
{
  printf("Command line:-%s",argv[0]);
  for(int i=1;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=3)
    {
      printf("Usage:- %s <q> <prec>.\n",argv[0]);
      return 0;
    }
  uint64_t q =atoi(argv[1]);
  if(q<3)
    {
      printf("q must be >=3.\n");
      return 0;
    }
  if((q%4)==2)
    {
      printf("There are no primitive characters of modulus %lu\n",q);
      return 0;
    }
      
  int64_t prec=atoi(argv[2]);
  if(prec<=0)
    {
      printf("Can't have non-positive precision.\n");
      return 0;
    }
  arb_t sqrt_q;
  dirichlet_group_t G;
  acb_dirichlet_hurwitz_precomp_t pre;
  acb_t tmp;
  acb_init(tmp);
  arb_set_d(acb_realref(tmp),0.5);
  //acb_dirichlet_hurwitz_precomp_init_num(pre,tmp,0,,prec);
  acb_clear(tmp);
  arb_init(sqrt_q);

  acb_t *w,*v;
  w=(acb_t *)malloc(sizeof(acb_t)*q);
  v=(acb_t *)malloc(sizeof(acb_t)*q);
  
  for(uint64_t i=0;i<q;i++)
    {
      acb_init(w[i]);
      acb_init(v[i]);
    }

  printf("w,v initialised.\n");
  

  dirichlet_char_t ch;

  arb_sqrt_ui(sqrt_q,q,prec);
  dirichlet_group_init(G,q);

  printf("Dirichlet Group initialised.\n");

  acb_t half,i_q;
  acb_init(half);acb_init(i_q);
  acb_set_d(half,0.5);
  for(uint64_t i=1;i<q;i++)
    if(co_prime(i,q))
      {
	acb_set_ui(i_q,i);
	acb_div_ui(i_q,i_q,q,prec);
	acb_dirichlet_hurwitz(v[i],half,i_q,prec);
      }
  acb_clear(half);
  acb_clear(i_q);
  
  acb_dirichlet_dft((acb_ptr) w,(acb_srcptr) v,G,prec);
  
  dirichlet_char_init(ch,G);
  dirichlet_char_first_primitive(ch,G);
  while(1==1)
    {
      uint64_t i=dirichlet_char_exp(G,ch);
      acb_div_arb(w[i],w[i],sqrt_q,prec);
      printf("%lu %lu: ",q,i);
      acb_printd(w[i],20);
      printf("\n");
      if(dirichlet_char_next_primitive(ch,G)<0)
	break;
    }
  dirichlet_char_clear(ch); 
  dirichlet_group_clear(G);


for(uint64_t i=0;i<q;i++)
  {
    acb_clear(w[i]);
    acb_clear(v[i]);
  }

arb_clear(sqrt_q);
//acb_dirichlet_hurwitz_precomp_clear(pre);  
  
  return 0;
}

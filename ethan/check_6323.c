#include "dirichlet.h"
#include "acb_dirichlet.h"


int main()
{
  dirichlet_group_t G;
  dirichlet_group_init(G,6323);


  dirichlet_char_t chi;
  dirichlet_char_init(chi,G);
  dirichlet_char_next_primitive(chi,G);

  arb_t ares;
  arb_init(ares);
  
  acb_t s,res;
  acb_init(s);acb_init(res);

  arb_set_d(acb_realref(s),0.5);
  arb_set_d(acb_imagref(s),0.0);

  do{
    acb_dirichlet_l(res,s,G,chi,100);
    acb_abs(ares,res,100);
    printf("%lu %lu ",dirichlet_index_char(G,chi),dirichlet_char_exp(G,chi));arb_printd(ares,30);printf("\n");
  } while (dirichlet_char_next_primitive(chi,G) >=0);
  
  return 0;
}

  

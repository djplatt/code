// print t and zeta'(1/2+t) for first zero above 10^9

#include "stdio.h"
#include "stdlib.h"
#include "flint/acb_dirichlet.h"

#define PREC (200)

int main()
{
  acb_t s;
  acb_struct *zetas;

  

  acb_init(s);
  zetas=(acb_struct *)malloc(sizeof(acb_t)*2);
  acb_init(zetas);
  acb_init(zetas+1);

  fmpz_t n;
  fmpz_init(n);
  fmpz_set_ui(n,2846548033ll);

  acb_dirichlet_zeta_zero(s,n,PREC);

  acb_dirichlet_zeta_jet(zetas,s,0,2,PREC);

  printf("s=");acb_printd(s,20);
  printf("\nzeta(s)=");acb_printd(zetas,20);
  printf("\nzeta'(s)=");acb_printd(zetas+1,20);
  acb_inv(s,zetas+1,PREC);
  printf("\n1/zeta'(s)=");acb_printd(s,20);printf("\n");

  return 0;
}

#include "stdlib.h"
#include "flint/acb_dirichlet.h"

int main()
{
  acb_t s;
  acb_struct *zeta;
  zeta=(acb_struct *)malloc(sizeof(acb_struct)*2);
  acb_init(s);
  acb_init(zeta);
  acb_init(zeta+1);
  arb_set_d(acb_realref(s),0.5);
  for(long unsigned int i=10000000000ll;i<10000000100ll;i++)
    {
      arb_set_ui(acb_imagref(s),i);
      acb_dirichlet_zeta_jet(zeta,s,0,2,135);
    }
  acb_printd(zeta+1,20);printf("\n");

  return 0;

}

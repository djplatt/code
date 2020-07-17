#include "malloc.h"
#include "stdlib.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../includes/mpfi_c.h"

#define debug printf("reached line %d\n",__LINE__);

mpfi_t **polyf;

void calc_polyf(unsigned int K)
{
  unsigned int k,coeff,k1_ptr=0,k2_ptr=2;
  if(!(polyf=(mpfi_t **) malloc(sizeof(mpfi_t *)*(K-1))))
    {
      printf("Failed to allocate memory for polyf. Exiting.\n");
      exit(1);
    }
  for(k=0;k<K-1;k++)
    {
      if(!(polyf[k]=(mpfi_t *) malloc(sizeof(mpfi_t)*(k+2))))
	{
	  printf("Failed to allocate memory for polyf. Exiting.\n");
	  exit(1);
	}
      for(coeff=0;coeff<k+2;coeff++)
	mpfi_init(polyf[k][coeff]);
    }

  mpfi_set_d(polyf[0][0],0.5);
  mpfi_set_d(polyf[0][1],1.0); // (f+1/2)
  for(k=1;k<K-1;k++)
    {
      mpfi_set_ui(polyf[k][k+1],1);
      mpfi_mul_d(polyf[k][0],polyf[k-1][0],0.5);
      for(coeff=1;coeff<=k;coeff++)
	{
	  mpfi_mul_d(polyf[k][coeff],polyf[k-1][coeff],coeff*2.0+0.5);
	  mpfi_add(polyf[k][coeff],polyf[k][coeff],polyf[k-1][coeff-1]);

	}

    }
}


int main()
{
  int i,j;
  mpfi_t *vec;
  mpfr_set_default_prec(100);
  /*
  vec=(mpfi_t *) malloc(sizeof(mpfi_t)*100);
  for(i=0;i<100;i++)
    {
      mpfi_init(vec[i]);
      mpfi_set_d(vec[i],i);
      mpfi_out_str(stdout,10,0,vec[i]);
    }
  */
  calc_polyf(10);
  for(i=0;i<9;i++)
    for(j=0;j<i+2;j++)
      mpfi_print(polyf[i][j]);
  exit(0);
}

#include "acb.h"
#include "inttypes.h"
#include "stdio.h"

int main(int argc, char **argv)
{
  FILE *infile=fopen("/projects/Zeros_of_Lfunctions/residues/residues_14.out","r");
if(!infile)
  {
    printf("residues_14.out not found. Exiting.\n");
    return 0;
  }
  arb_t tmp,rho_bit,quarter,term,total;
  arb_init(tmp);arb_init(rho_bit);arb_init(quarter);
  arb_init(term);arb_init(total);
  arb_set_d(quarter,0.25);
  acb_t residue;
  acb_init(residue);
  uint64_t prec=200;

  uint64_t i;
  for(i=0;i<2702;i++) // gamma < exp(exp(2))
    {
      arb_load_file(tmp,infile);
      printf("%4d gamma = ",i+1);arb_printd(tmp,30);printf("\n");
      arb_sqr(rho_bit,tmp,prec);
      arb_add(tmp,rho_bit,quarter,prec);
      arb_sqrt(rho_bit,tmp,prec); // sqrt(1/4+gamma^2)
      arb_load_file(acb_realref(residue),infile);
      arb_load_file(acb_imagref(residue),infile);
      printf("%4d 1/zeta'(1/2+i gamma) = ",i+1);acb_printd(residue,15);printf("\n");
      acb_abs(tmp,residue,prec); // |1/zeta'(1/2+rho)|
      arb_div(term,tmp,rho_bit,prec); // 1/(sqrt(1/4+gamma^2)|zeta'(1/2+igamma)|)
      printf("%4d term = ",i+1);arb_printd(term,30);printf("\n");
      arb_add(total,total,term,prec);
    }
  printf("Sum = ");arb_printd(total,30);printf("\n");
  return 0;
}


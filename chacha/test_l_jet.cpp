#include "inttypes.h"
#include "dirichlet.h"
#include "acb_dirichlet.h"

#define OP_WIDTH (10)

int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("Usage:- %s <log 2 epsilon>\n",argv[0]);
      return 0;
    }

  int64_t l2delta=atol(argv[1]);
  acb_t z,ll,z1;arb_t dt,tmp;
  acb_init(ll);
  acb_init(z);acb_init(z1);
  arb_set_d(acb_realref(z),0.5);
  arb_set_d(acb_imagref(z),1.5);
  arb_init(dt);
  arb_set_ui(dt,1);
  arb_mul_2exp_si(dt,dt,l2delta);
  printf("epsilon = ");arb_printd(dt,OP_WIDTH);printf("\n");
  acb_set(z1,z);
  arb_add_error(acb_realref(z1),dt);
  arb_add_error(acb_imagref(z1),dt);
  arb_init(tmp);

  dirichlet_group_t dg;
  dirichlet_char_t chi;
  uint64_t q=123;

  dirichlet_group_init(dg,q);
  dirichlet_char_init(chi,dg);
  dirichlet_char_first_primitive(chi,dg);

  acb_ptr res1,res2;
  res1=(acb_ptr)malloc(sizeof(acb_t)*4);
  res2=(acb_ptr)malloc(sizeof(acb_t)*4);

  for(uint64_t n=0;n<4;n++)
    {acb_init(res1+n);acb_init(res2+n);}

  printf("exact\n");
  for(uint64_t n=0;n<4;n++)
    {
      printf("\\begin{table}[h]\n\\centering\n\\begin{tabular}{c c}\nPrecision (bits) & $L^{(%lu)}(z)$\\\\\n\\hline\n",n);
      for(int64_t prec=4;prec<=128;prec+=prec)
	{
	  acb_dirichlet_l_jet(res1,z,dg,chi,0,n+1,prec);
	  acb_mul_ui(res1+2,res1+2,2,prec);
	  acb_mul_ui(res1+3,res1+3,6,prec);
	  printf("$%ld$ & $",prec);
	  acb_printd(res1+n,OP_WIDTH);printf("$ \\\\\n");
	}
      printf("\\end{tabular}\n\\caption{With z exact.}\n\\end{table}\n");
    }

  printf("inexact\n");
  for(uint64_t n=0;n<4;n++)
    {
      printf("\\begin{table}[h]\n\\centering\n\\begin{tabular}{c c}\nPrecision (bits) & $L^{(%lu)}(z)$\\\\\n\\hline\n",n);
      for(int64_t prec=4;prec<=128;prec+=prec)
	{
	  acb_dirichlet_l_jet(res1,z1,dg,chi,0,n+1,prec);
	  acb_mul_ui(res1+2,res1+2,2,prec);
	  acb_mul_ui(res1+3,res1+3,6,prec);
	  printf("$%ld$ & $",prec);
	  acb_printd(res1+n,OP_WIDTH);printf("$ \\\\\n");
	}
      printf("\\end{tabular}\n\\caption{With $z\\pm\\delta$.}\n\\end{table}\n");
    }

  /*


  acb_dirichlet_l_jet(res2,z1,dg,chi,0,4,prec);
      acb_mul_ui(res2+2,res2+2,2,prec);

  acb_mul_ui(res2+3,res2+3,6,prec);



  printf("\ninexact\n");
  for(uint64_t n=0;n<4;n++)
    {acb_printd(res2+n,OP_WIDTH);printf("\n");}
  acb_div(ll,res2+2,res2+1,prec);
  printf("L''/L' = ");acb_printd(ll,OP_WIDTH);

  printf("\n\ncomputed inexeact\n");

  acb_set(res1+3,res2+3);
  for(uint64_t n=3;n>0;n--)
    {
      acb_abs(tmp,res1+n,prec);
      arb_mul(tmp,tmp,dt,prec);
      arb_add_error(acb_realref(res1+n-1),tmp);
      arb_add_error(acb_imagref(res1+n-1),tmp);
    }
  for(uint64_t n=0;n<4;n++)
    {acb_printd(res1+n,OP_WIDTH);printf("\n");}

  acb_div(ll,res1+2,res1+1,prec);
  printf("L''/L' = ");acb_printd(ll,OP_WIDTH);printf("\n");
    }
  */
  return 0;
}

/*

compute int_X^2X |psi(x)-x|2 dx/X^2
low memory footprint version
using ARB

*/
  
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "pari/pari.h"
#include "inttypes.h"
#include "arb.h"


void term(arb_ptr res, uint64_t a, arb_ptr psa, int64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_sub_ui(tmp,psa,2*a+1,prec);
  arb_mul(res,tmp,psa,prec);
}
//double res=psa*(psa-2*a-1);
//return res;

// return p if n=p^k, 0 otherwise
uint64_t myisprimepower(uint64_t n)
{
  pari_sp av;
  av=avma;
  GEN nn,p[1];
  nn=uutoi(0,n); // nn=n
  p[0]=uutoi(0,0); // p[0]=0
  isprimepower(nn,p);
  uint64_t res=itos(p[0]);
  avma=av;
  //printf("isprimepower(%lu) returning %lu\n",n,res);
  return res;
}

void out(arb_t I, uint64_t X, int64_t prec)
{
  static bool init=false;
  static arb_t tmp,tmp1,tmp2;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
    }
  arb_set_ui(tmp,X*7);
  arb_div_ui(tmp1,tmp,3,prec);
  arb_div_ui(tmp,I,X,prec);
  arb_div_ui(tmp2,tmp,X,prec);
  arb_add(tmp,tmp1,tmp2,prec);
  arb_printd(tmp,30);
}

int main(int argc, char **argv)
{

  if(argc!=4)
    {
      printf("Usage:- %s <Max X> <step> <prec>.\n",argv[0]);
      return 0;
    }
  uint64_t MAX_X=atol(argv[1]);
  uint64_t step=atol(argv[2]);
  int64_t prec=atol(argv[3]);

  pari_init(1<<28,1000);

  uint64_t X=1;
  arb_t I,PsiX,Psi2X1,tmp;
  arb_init(I);arb_init(PsiX);arb_init(Psi2X1);arb_init(tmp);
  term(I,1,PsiX,prec);
  //double I=term(1,0.0); // I now contains int_1^2 psi(t)^2 -2 t psi(t) dt 
  //double PsiX=0.0,Psi2X1=0.0;
  //printf("%lu ",X);out(I,X,prec);


  // on entry X=1, PsiX=psi(0), Psi2X1=psi(1) 
  while(X<MAX_X)
    {
      X++;
      term(tmp,X-1,PsiX,prec);
      arb_sub(I,I,tmp,prec);
      uint64_t p=myisprimepower(X);
      if(p) 
	{
	  arb_log_ui(tmp,p,prec);
	  arb_add(PsiX,PsiX,tmp,prec);
	}
      p=myisprimepower(X+X-2);
      if(p) 
	{
	  arb_log_ui(tmp,p,prec);
	  arb_add(Psi2X1,Psi2X1,tmp,prec);
	}
      term(tmp,X+X-2,Psi2X1,prec);
      arb_add(I,I,tmp,prec);
      p=myisprimepower(X+X-1);
      if(p) 
	{
	  arb_log_ui(tmp,p,prec);
	  arb_add(Psi2X1,Psi2X1,tmp,prec);
	}
      term(tmp,X+X-1,Psi2X1,prec);
      arb_add(I,I,tmp,prec);
      if((X%step)==0)
	{
	  printf("%lu ",X);out(I,X,prec);printf(" ");arb_printd(Psi2X1,20);printf("\n");
	  fflush(stdout);
	}
    }

  return 0;
}

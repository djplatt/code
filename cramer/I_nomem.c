/*

compute int_X^2X |psi(x)-x|2 dx/X^2
low memory footprint version

*/
  
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "pari/pari.h"

#define MAX_X (100000000)

double term(uint64_t a, double psa)
{
  double res=psa*(psa-2*a-1);
  //printf("term(%lu,%e) returning %e\n",a,psa,res);
  return res;
}

// return p if n=p^k, 0 otherwise
uint64_t myisprimepower(uint64_t n)
{
  pari_sp av;
  av=avma;
  GEN nn,p[1];
  nn=uutoi(n>>32,n&0xffffffff); // nn=n
  p[0]=uutoi(0,0); // p[0]=0
  isprimepower(nn,p);
  uint64_t res=itos(p[0]);
  avma=av;
  //printf("isprimepower(%lu) returning %lu\n",n,res);
  return res;
}

int main(int argc, char **argv)
{

  pari_init(1<<28,1000);

  uint64_t X=1;
  double I=term(1,0.0); // I now contains int_1^2 psi(t)^2 -2 t psi(t) dt 
  double PsiX=0.0,Psi2X1=0.0;
  printf("%lu %e\n",X,7.0*X/3.0+I/(X*X));

  // on entry X=1, PsiX=psi(0), Psi2X1=psi(1) 
  while(X<MAX_X)
    {
      X++;
      I-=term(X-1,PsiX); 
      uint64_t p=myisprimepower(X);
      if(p) PsiX+=log(p);
      p=myisprimepower(X+X-2);
      if(p) Psi2X1+=log(p); // Psi2X1=psi(2X-2)
      I+=term(X+X-2,Psi2X1); // add int 2X-2..2X-1
      p=myisprimepower(X+X-1);
      if(p) Psi2X1+=log(p); // Psi2X1=psi(2X-1)
      I+=term(X+X-1,Psi2X1); // add int 2X-1..2X
      double X2=X;
      X2*=X;
      printf("%lu %e\n",X,7.0*X/3.0+I/X2);
    }

  return 0;
}


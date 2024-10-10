#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct{
 double k;
 double l;
 double k_l;} k_l_t;

void print_k_l(k_l_t k_l)
{
  printf("(k,l)=(%16.14f,%16.14f) : k+l-1/2 = %16.14f.\n",k_l.k,k_l.l,k_l.k_l);
}

k_l_t B(k_l_t inp)
{
  k_l_t op;
  op.k=inp.l-0.5;
  op.l=inp.k+0.5;
  op.k_l=inp.k_l;
  return op;
}

k_l_t A(k_l_t inp)
{
  k_l_t op;
  double num=2.0*(1.0+inp.k);
  op.k=inp.k/num;
  op.l=0.5+inp.l/num;
  op.k_l=op.k+op.l-0.5;
  return op;
}

k_l_t AB(k_l_t inp)
{
  return A(B(inp));
}

k_l_t ABA(k_l_t inp)
{
  return AB(A(inp));
}

k_l_t ABAA(k_l_t inp)
{
  return ABA(A(inp));
}

int main()
{
  k_l_t init;
  init.k=0.0;
  init.l=55.0/84.0;
  init.k_l=init.k+init.l-0.5;
  init=B(init);
  printf("starting with: ");print_k_l(init);

  for(uint64_t n=0;n<11;n++)
    {
      print_k_l(init);
      k_l_t ab=AB(init);printf("AB -> ");print_k_l(ab);
      k_l_t aba=ABA(init);printf("ABA -> ");print_k_l(aba);
      k_l_t abaa=ABAA(init);printf("ABAA -> ");print_k_l(abaa);
      if(ab.k_l<aba.k_l)
	{
	  if(ab.k_l<abaa.k_l)
	    {printf("AB\n");init=ab;}
	  else
	    {printf("ABAA\n");init=abaa;}
	}
      else
	{
	  if(aba.k_l<abaa.k_l)
	    {printf("ABA\n");init=aba;}
	  else
	    {printf("ABAA\n");init=abaa;}
	}
    }
  
  print_k_l(init);
  
  return 0;
}
  

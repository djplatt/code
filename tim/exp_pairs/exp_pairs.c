#include <inttypes.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <flint/arb.h>

arb_t half;

void A(arb_t new_k, arb_t new_l, arb_t k, arb_t l, int prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
    }
  arb_add_ui(tmp1,k,1,prec);
  arb_mul_2exp_si(tmp1,tmp1,1); // 2(1+k)
  arb_div(new_k,k,tmp1,prec);
  arb_div(tmp2,l,tmp1,prec);
  arb_add(new_l,tmp2,half,prec);
}

void B(arb_t new_k, arb_t new_l, arb_t k, arb_t l, int prec)
{
  static bool init=false;
  static arb_t old_k;
  if(!init)
    {
      init=true;
      arb_init(old_k);
    }
  arb_set(old_k,k);
  arb_sub(new_k,l,half,prec);
  arb_add(new_l,old_k,half,prec);
}

void ABA0(arb_t new_k, arb_t new_l, arb_t k, arb_t l, int prec)
{
  B(new_k,new_l,k,l,prec);
  A(new_k,new_l,new_k,new_l,prec);
}

void ABA1(arb_t new_k, arb_t new_l, arb_t k, arb_t l, int prec)
{
  A(new_k,new_l,k,l,prec);
  B(new_k,new_l,new_k,new_l,prec);
  A(new_k,new_l,new_k,new_l,prec);
}

void ABA2(arb_t new_k, arb_t new_l, arb_t k, arb_t l, int prec)
{
  A(new_k,new_l,k,l,prec);
  A(new_k,new_l,new_k,new_l,prec);
  B(new_k,new_l,new_k,new_l,prec);
  A(new_k,new_l,new_k,new_l,prec);
}

bool smaller_p (arb_t a, arb_t b, int64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_sub(tmp,a,b,prec);
  return arb_is_negative(tmp);
}
			  

void greedy(arb_t best_k, arb_t best_l, arb_t k, arb_t l, int64_t prec)
{
  static bool init=false;
  static arb_t ks[3],ls[3],k_ls[3];
  if(!init)
    {
      init=true;
      for(uint64_t i=0;i<3;i++)
	{
	  arb_init(ks[i]);
	  arb_init(ls[i]);
	  arb_init(k_ls[i]);
	}
    }
  ABA0(ks[0],ls[0],k,l,prec);
  ABA1(ks[1],ls[1],k,l,prec);
  ABA2(ks[2],ls[2],k,l,prec);
  for(uint64_t i=0;i<3;i++)
    arb_add(k_ls[i],ks[i],ls[i],prec);
  if(smaller_p(k_ls[1],k_ls[0],prec)) // ABA beat AB
    {
      if(smaller_p(k_ls[2],k_ls[1],prec)) // ABAA beat ABA
	{
	  printf("Applying ABAA.\n");
	  arb_set(best_k,ks[2]);
	  arb_set(best_l,ls[2]);
	  return;
	}
      else
	{
	  printf("Applying ABA.\n");
	  arb_set(best_k,ks[1]);
	  arb_set(best_l,ls[1]);
	  return;
	}
    }
  else
    {
      if(smaller_p(k_ls[2],k_ls[0],prec))
	{
	  printf("Applying ABAA.\n");
	  arb_set(best_k,ks[2]);
	  arb_set(best_l,ls[2]);
	  return;
	}
      else
	{
	  printf("Applying AB.\n");
	  arb_set(best_k,ks[0]);
	  arb_set(best_l,ls[0]);
	  return;
	}
    }
}

int main()
{

  int64_t prec=1000;
  arb_init(half);
  arb_set_d(half,0.5);

  arb_t k0,l0,k,l;
  arb_init(k0);arb_init(l0);arb_init(k);arb_init(l);

  arb_set_d(k,0.0);
  arb_set_d(l,1.0);

  
  
  for(uint64_t n=0;n<100;n++)
    greedy(k,l,k,l,prec);      

  //ABA0(k,l,k,l,prec);
  //ABA1(k,l,k,l,prec);
  //ABA1(k,l,k,l,prec);
  //ABA1(k,l,k,l,prec);
  //ABA2(k,l,k,l,prec); // ABA^3BA^2BA^2B = Phillips

  arb_add(k,k,l,prec);
  arb_sub(k,k,half,prec);

  printf("k+l-1/2 = ");arb_printd(k,20);printf("\n");
  

  return 0;
}

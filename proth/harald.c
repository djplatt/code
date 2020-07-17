#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <math.h>
#include <gmp.h>

#define Nsp 16
int smallprime[Nsp] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53};
mpz_t smallprimebig[Nsp];

int jac(int a, int b)
/* computes Jacobi symbol (a/b) */
/* preconditions: a,b>0 */
{
  int r,s;
  int prod;

  if(b==1)
    return 1;
  a = a%b;
  if(a==0)
    return 0;
  
  if(!(b%2)) {
    if(!(a%2))
      return 0;  
    while(!(b%2))
      b/=2;
  }
  
  s=0; 
  while(!(a%2)) {
    a/=2; s++;
  }
  if(s%2) {
    if(((b%8)==1) || ((b%8)==7))
      prod = 1;
    else
      prod = -1;
  } else prod = 1;

  /* At this stage, prod = (-1)^((1/8) (b^2-1)) if v_2(a) was odd,
     v_2(b) = 0 and b\equiv +- 1 mod 8; otherwise prod = 1.
     Both a and b have been divided by the powers of 2 they contained. */

  return prod*((((a-1)%4) && ((b-1)%4)) ? -1 : 1)*jac(b,a);
  /*by quadratic reciprocity */
}
                        
int jacill(int a, mpz_t b)
/* computes Jacobi symbol (a/b) */
/* preconditions: a,b>0 */
/* accepts int a and mpz_t b as inputs; calls jac */
{
  unsigned long int r, bmoda;
  int s;
  int prod;
  mpz_t scratch;
  int ojac;
  
  if(!mpz_cmp_si(b,1))
    return 1;
  if(a==0)
    return 0;

  mpz_init(scratch);  
  if(mpz_divisible_ui_p(b,2)) {
    if(!(a%2)) {
      return 0;  
    }
    while(mpz_divisible_ui_p(b,2)) {
      mpz_divexact_ui(scratch,b,2);
      mpz_set(b,scratch);
    }
  }
   
  s=0; 
  while(!(a%2)) {
    a/=2; s++;
  }
  if(s%2) {
    r = mpz_mod_ui(scratch,b,8);
    if(((r%8)==1) || ((r%8)==7))
      prod = 1;
    else
      prod = -1;
  } else prod = 1;
  
  /* At this stage, prod = (-1)^((1/8) (b^2-1)) if v_2(a) was odd,
     v_2(b) = 0 and b\equiv +- 1 mod 8; otherwise prod = 1.
     Both a and b have been divided by the powers of 2 they contained. */

  r = mpz_mod_ui(scratch,b,4);
  bmoda = mpz_mod_ui(scratch,b,a);
  mpz_clear(scratch);
  return prod*((((a-1)%4) && ((r-1)%4)) ? -1 : 1)*jac(bmoda,a);
  /*by quadratic reciprocity */
}

int isproth(mpz_t p)
/* Precondition: p of the form k2^n+1, k odd, k<2^n */
/* Returns 1 if p is guaranteed to be a prime */        
/* Returns 0 otherwise */
{
  mpz_t res,halfpminusone,pminusone;
  int i;

  for(i=0; i<Nsp; i++)
    if(jacill(smallprime[i],p)==-1)
      break;
  if(i==Nsp)
    return 0;
  mpz_init(res); mpz_init(halfpminusone); mpz_init(pminusone);
  mpz_fdiv_q_ui(halfpminusone,p,2);
  mpz_powm(res,smallprimebig[i],halfpminusone,p);
  mpz_clear(halfpminusone); 
  mpz_sub_ui(pminusone,p,1);
  if(!mpz_cmp(res,pminusone)) {
    mpz_clear(res); mpz_clear(pminusone);
    return 1;
  } else {
    mpz_clear(res); mpz_clear(pminusone);
    return 0;
  }
}

main(int argc, char *argv[])
{
  mpz_t p,twoton,twicethat,twototwon,scratch,step,os,newp;
  int i;
  int n;
  

  if(argc>1)
    n=atoi(argv[1]);
  else n=30;
  
  mpz_init(os);
  if(argc>2)
    mpz_set_str(os,argv[2],10);
  else mpz_set_ui(os,4000000000000000000);
  /*    os = 4e18;*/
  /*  mpz_out_str(stdout,10,os);*/
  for(i=0; i<Nsp; i++)
    mpz_init_set_si(smallprimebig[i],smallprime[i]);
 
  mpz_init(p); mpz_init(twoton); mpz_init(twicethat); mpz_init(scratch);
  mpz_init(twototwon); mpz_init(step); mpz_init(newp);
  mpz_ui_pow_ui(twoton,2,n);
  mpz_ui_pow_ui(twicethat,2,n+1);
  mpz_ui_pow_ui(twototwon,2,2*n);
  mpz_fdiv_r(scratch,os,twicethat);
  mpz_sub(step,os,scratch);
  mpz_add_ui(p,twoton,1);
  while(mpz_cmp(p,twototwon)<0) { 
    if(isproth(p)) {
      mpz_add(newp,p,step);
      //mpz_out_str(stdout,10,p); printf("\n");
      break;
    } else {
      mpz_add(scratch,p,twicethat);
      mpz_set(p,scratch);
    }
  }
  
  while(mpz_cmp(newp,twototwon)<0) {
    while(!isproth(newp) && mpz_cmp(newp,p)>0) {
      mpz_sub(scratch,newp,twicethat);
      mpz_set(newp,scratch);
    }
    if(mpz_cmp(newp,p)<=0) {
      fprintf(stderr,"I kill myself at ");
      mpz_out_str(stderr,10,p); fprintf(stderr,"\n");
      exit(1);
    } 
    //mpz_out_str(stdout,10,newp); printf("\n");
    mpz_set(p, newp); mpz_add(newp, p, step);  
  }
}

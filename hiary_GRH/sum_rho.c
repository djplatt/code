#include "stdlib.h"
#include "flint/acb_dirichlet.h"
#include "inttypes.h"
#include "stdbool.h"
#define OP_ACC ((int64_t) 101)
#define ZERO_LEN ((uint64_t) 13)


void XGCD(int64_t *d, 
	  int64_t *s, 
	  int64_t *t, 
	  int64_t a, 
	  int64_t b)
{
   int64_t  u, v, u0, v0, u1, v1, u2, v2, q, r;

   int64_t aneg = 0, bneg = 0;

   if (a < 0) {
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      b = -b;
      bneg = 1;
   }

   u1=1; v1=0;
   u2=0; v2=1;
   u = a; v = b;

   while (v != 0) {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 =  u1 - q*u2;
      v2 = v1- q*v2;
      u1 = u0;
      v1 = v0;
   }

   if (aneg)
      u1 = -u1;

   if (bneg)
      v1 = -v1;

   d[0] = u;
   s[0] = u1;
   t[0] = v1;
}
   
// taken from NTL:
int64_t InvMod(int64_t a, int64_t n)
{
   int64_t d, s, t;

   XGCD(&d, &s, &t, a, n);
   if (d != 1) return -1;
   if (s < 0)
      return s + n;
   else
      return s;
}


uint64_t gcd (uint64_t a, uint64_t b)
/* Euclid algorithm gcd */
{
  uint64_t c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
};

int co_prime(uint64_t a, uint64_t b)
{
  return(gcd(a,b)==1);
};

// exp(2 pi i x)
void arb_e(acb_t res, arb_t x, uint64_t prec)
{
  static arb_t two_pi;
  static acb_t temp;
  static bool init=false;

  if(!init)
    {
      init=true;
      arb_init(two_pi);
      arb_const_pi(two_pi,prec);
      arb_mul_2exp_si(two_pi,two_pi,1);
      acb_init(temp);
    }

  arb_mul(acb_imagref(temp),two_pi,x,prec);
  acb_exp(res,temp,prec);
}
  
void in_bytes(arb_ptr t, FILE *infile, uint64_t prec)
{
  uint64_t a;
  uint32_t b;
  uint8_t c;

  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (a). Exiting.\n");
      exit(0);
    }
  //printf("a=%lu\n",a);
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (b). Exiting.\n");
      exit(0);
    }
  //printf("b=%u\n",b);
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (c). Exiting.\n");
      exit(0);
    }
  //printf("c=%u\n",c);
  arb_set_ui(t,c);
  arb_mul_2exp_si(t,t,32);
  arb_add_ui(t,t,b,prec);
  arb_mul_2exp_si(t,t,64);
  arb_add_ui(t,t,a,prec);
  arb_mul_2exp_si(t,t,-OP_ACC);
}

bool do_rho1(arb_t res, arb_t full_sum, arb_t gam, int64_t prec)
{
  static bool init=false;
  static arb_t pm1,qtr,rho,rho2;
  if(!init)
    {
      init=true;
      arb_init(pm1);
      arb_set_ui(pm1,1);
      arb_mul_2exp_si(pm1,pm1,-OP_ACC-1);
      arb_init(qtr);
      arb_set_d(qtr,0.25);
      arb_init(rho);
      arb_init(rho2);
    }

  arb_add(rho,gam,pm1,prec);
  arb_sqr(rho2,rho,prec);
  arb_add(rho,rho2,qtr,prec);
  arb_inv(rho2,rho,prec);
  arb_add(res,res,rho2,prec);
  arb_add(rho,res,rho2,prec);
  arb_sub(rho2,rho,full_sum,prec);
  return arb_is_positive(rho2);
}

void do_rho(arb_t res, arb_t gam, int64_t prec)
{
  static bool init=false;
  static arb_t pm1,qtr,rho,rho2;
  if(!init)
    {
      init=true;
      arb_init(pm1);
      arb_set_ui(pm1,1);
      arb_mul_2exp_si(pm1,pm1,-OP_ACC-1);
      arb_init(qtr);
      arb_set_d(qtr,0.25);
      arb_init(rho);
      arb_init(rho2);
    }

  arb_add(gam,gam,pm1,prec);
  arb_sqr(rho2,gam,prec);
  arb_add(rho,rho2,qtr,prec);
  arb_inv(res,rho,prec);
}

int main(int argc, char** argv)
{
  printf("Command line:- ");
  for(uint64_t i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=3)
    {
      printf("Fatal error in main: usage %s <prec> <zeros file>. Exiting.\n",argv[0]);
      exit(0);
    }
  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Fatal error in main: failed to open zeros file for binary input. Exiting.\n");
      exit(0);
    }
  uint64_t prec=atol(argv[1]);
  
  uint64_t q,index,num_zeros;
  if(fread(&q,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading q from %s. Exiting.\n",argv[2]);
      exit(0);
    }

  
  arb_t t,del_t;
  arb_init(t);
  arb_init(del_t);
  arb_t sum,sum1;arb_init(sum);arb_init(sum1);

  while(fread(&index,sizeof(uint64_t),1,infile)==1)
    {
      if(fread(&num_zeros,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Error reading num_zeros from %s. Exiting.\n",argv[2]);
	  exit(0);
	}
      if(index!=InvMod(index,q)) // not real
	{
	  if(fseek(infile,num_zeros*ZERO_LEN,SEEK_CUR)!=0)
	    {
	      printf("Fatal error skipping zeros %lu %lu. Exiting.\n",q,index);
	      exit(0);
	    }
	  continue;
	}
      arb_t *zeros=(arb_t *)malloc(sizeof(arb_t)*num_zeros);
      arb_t *inv_rhos=(arb_t *)malloc(sizeof(arb_t)*num_zeros);
      arb_t *psums=(arb_t *)malloc(sizeof(arb_t)*num_zeros);
      if(q==3)
	arb_set_ui(t,8);
      else
	arb_set_ui(t,0);
      in_bytes(del_t,infile,prec);
      arb_init(zeros[0]);
      arb_add(zeros[0],t,del_t,prec);
      arb_init(inv_rhos[0]);
      do_rho(inv_rhos[0],zeros[0],prec);
      arb_init(psums[0]);
      arb_set(psums[0],inv_rhos[0]);
      for(uint64_t z=1;z<num_zeros;z++)
	{
	  in_bytes(del_t,infile,prec);
	  arb_add(t,t,del_t,prec);
	  arb_init(zeros[z]);
	  arb_set(zeros[z],t);
	  arb_init(inv_rhos[z]);
	  do_rho(inv_rhos[z],zeros[z],prec);
	  arb_init(psums[z]);
	  arb_add(psums[z],psums[z-1],inv_rhos[z],prec);
	}
      //printf("Zero %lu at ",num_zeros);arb_printd(zeros[num_zeros-1],20);printf("\n");
      for(uint64_t z=0;z<num_zeros;z++)
	{
	  printf("Zero %lu at 1/2+ i ",z+1);arb_printd(zeros[z],10);
	  printf(" contributes ");arb_printd(inv_rhos[z],10);
	  printf(" sums to ");arb_printd(psums[z],10);printf("\n");
	}
      for(uint64_t z=0;z<num_zeros;z++)
	{
	  arb_clear(zeros[z]);
	  arb_clear(inv_rhos[z]);
	  arb_clear(psums[z]);
	}
      free(zeros);
      free(inv_rhos);
      free(psums);
    }
  return 0;
}

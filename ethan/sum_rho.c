/*

File: sum_gamma1.c

Created: 9 Feb 2023

Version: 1.0

Last Modified: 

Dialect: C

Requires: ARB


V 1.0 Initial implementation


Build instructions: gcc -O2 -larb

By: DJ Platt
    Bristol University

Copyright 2023.*/

#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "acb.h"
#include "stdbool.h"

#define OP_ACC ((int64_t) 101)


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


void in_bytes(arb_ptr t, FILE *infile, uint64_t prec)
{
  uint64_t a;
  uint32_t b;
  uint8_t c;
  int res;

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
  arb_set_ui(t,c);
  arb_mul_2exp_si(t,t,32);
  arb_add_ui(t,t,b,prec);
  arb_mul_2exp_si(t,t,64);
  arb_add_ui(t,t,a,prec);
  arb_mul_2exp_si(t,t,-OP_ACC);
}

void do_rho(arb_ptr res, arb_ptr rho, uint64_t prec)
{
  static bool init=false;
  static arb_t quarter;
  if(!init)
    {
      init=true;
      arb_init(quarter);
      arb_set_d(quarter,0.25);
    }
  arb_mul(res,rho,rho,prec);
  arb_add(res,res,quarter,prec);
  arb_sqrt(res,res,prec);
  arb_inv(res,res,prec);
}
 
int main(int argc, char **argv)
{
  int i;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
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
  
  uint64_t q,index,num_zeros,index1,num_zeros1,num_prims=0;
  if(fread(&q,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading q from %s. Exiting.\n",argv[2]);
      exit(0);
    }

  arb_t tot,res,rho,t,del_t;
  arb_init(tot);arb_init(res);arb_init(rho);arb_init(t);arb_init(del_t);



  while(fread(&index,sizeof(uint64_t),1,infile)==1)
    {
      num_prims++;
      //printf("Doing index %lu\n",index);
      if(fread(&num_zeros,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Error reading num_zeros from %s. Exiting.\n",argv[2]);
	  exit(0);
	}
      //printf("Doing %lu zeros\n",num_zeros);
      arb_set_ui(tot,0);
      uint64_t z;
      if(q==3)
	arb_set_ui(t,8);
      else
	arb_set_ui(t,0);
      for(z=0;z<num_zeros;z++)
	{
	  in_bytes(del_t,infile,prec);
	  arb_add(t,t,del_t,prec); // exact
	  arb_set(rho,t);
	  arb_add_error_2exp_si(rho,-OP_ACC-1);
	  do_rho(res,rho,prec); // returns 1/rho
	  if((z==0)&&(index==3837))
	    {
	      arb_printd(rho,20);
	      printf("\n");
	    }
	      
	  arb_add(tot,tot,res,prec);
	}
      uint64_t conj_index=InvMod(index,q);
      if(conj_index==index) // a real character
	{
	  arb_mul_2exp_si(tot,tot,1);
	  printf("q %lu index %lu %lu sum (1/4+gamma^2)^(-1/2) = ",q,index,index);arb_printd(tot,20);printf("\n");
	  continue;
	}
      num_prims++;
      if(fread(&index1,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Error reading index1 from %s. Exiting.\n",argv[2]);
	  exit(0);
	}
      //printf("Conjugate index read was %lu\n",index1);
      if(index1!=conj_index)
	{
	  printf("Conj index read (%lu) does not match expected (%lu). Exiting.\n",index1,conj_index);
	  exit(0);
	}
      if(fread(&num_zeros1,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Error reading num_zeros from %s. Exiting.\n",argv[2]);
	  exit(0);
	}
      //printf("Doing %lu conjugate zeros.\n",num_zeros1);
      arb_set_ui(t,0);
      for(z=0;z<num_zeros1;z++)
	{
	  in_bytes(del_t,infile,prec);
	  arb_add(t,t,del_t,prec); // exact
	  arb_set(rho,t);
	  arb_add_error_2exp_si(rho,-OP_ACC-1);
	  do_rho(res,rho,prec); // returns 1/rho
	  if((z==0)&&(index==3837))
	    {
	      arb_printd(rho,20);
	      printf("\n");
	    }
	  arb_add(tot,tot,res,prec);
	}
      printf("q %lu index %lu %lu sum (1/4+gamma^2)^(-1/2) = ",q,index,conj_index);arb_printd(tot,20);printf("\n");
      
    }

  return(0);
}

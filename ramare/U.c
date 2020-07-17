/*

File: U.c

Created: 21 January 2014

Version: 1.0

Last Modified: 

Dialect: C

Requires: GMP
          MPFR
          MPFI

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

V 1.0 Initial implementation


Build instructions: gcc -O2 -lmpfi

By: DJ Platt
    Bristol University

Copyright 2012.

The author is funded by the UK ESPRC. */

#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "inttypes.h"
#include "../includes/mpfi_c.h"

#define OP_ACC (101)
#include "inttypes.h"


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



mpfi_t pm1;
void init_in_bytes()
{
  mpfi_t tmp;
  mpfi_init(tmp);
  mpfi_init(pm1);
  mpfi_set_ui(pm1,1);
  mpfi_neg(tmp,pm1);
  mpfi_put(pm1,tmp);
  mpfi_div_2ui(pm1,pm1,OP_ACC+1);
  mpfi_clear(tmp);
}

void in_bytes(mpfi_ptr t, FILE *infile)
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
  mpfi_set_ui(t,c);
  mpfi_mul_2ui(t,t,32);
  mpfi_add_ui(t,t,b);
  mpfi_mul_2ui(t,t,64);
  mpfi_add_ui(t,t,a);
  mpfi_div_2ui(t,t,OP_ACC);
}


// read the next imaginary part of rho into del_t
inline void next_rho(mpfi_ptr del_t, FILE *infile)
{
  in_bytes(del_t,infile);
}

mpfi_t rho_tmp;
void do_rho(mpfi_ptr res, mpfi_ptr rho)
{
  mpfi_sqr(res,rho);
  mpfi_add_d(rho_tmp,res,0.25);
  mpfi_inv(res,rho_tmp);
}

mpfi_t del_t,t,rho,res;
void do_zeros(mpfi_ptr tot, uint64_t num_zeros, uint64_t q, FILE *infile)
{
  uint64_t z;
  if(q==3)
    mpfi_set_ui(t,8);
  else
    mpfi_set_ui(t,0);
  for(z=0;z<num_zeros;z++)
    {
      next_rho(del_t,infile);
      mpfi_add(t,t,del_t); // exact
      mpfi_add(rho,t,pm1); // +/- OP_ACC
      do_rho(res,rho);
      //mpfi_print_str("do_rho returned ",res);
      mpfi_add(tot,tot,res);
    }
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
      printf("Fatal error in main: usage %s <prec> <zeros file>. Exiting.\n");
      exit(0);
    }
  mpfr_set_default_prec(atoi(argv[1]));
  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Fatal error in main: failed to open zeros file for binary input. Exiting.\n");
      exit(0);
    }

  uint64_t q,index,num_zeros,index1,num_zeros1,num_prims=0;
  if(fread(&q,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading q from %s. Exiting.\n",argv[2]);
      exit(0);
    }

  mpfi_t tot,tot1;
  mpfi_init(res);mpfi_init(tot);mpfi_init(del_t);mpfi_init(t);
  mpfi_init(rho);mpfi_init(tot1);
  mpfi_set_ui(tot1,0);

  mpfi_init(rho_tmp);

  init_in_bytes();

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
      mpfi_set_ui(tot,0);
      do_zeros(tot,num_zeros,q,infile);
      uint64_t conj_index=InvMod(index,q);
      if(conj_index==index) // a real character
	{
	  // We have sum for Im rho in [0,H]
	  // double it to get |Im rho| <=H
	  mpfi_add(tot1,tot1,tot);
	  mpfi_add(tot1,tot1,tot);
	  //printf("%lu %lu %lu",q,index,num_zeros*2);
	  //mpfi_print_str(" sum ",tot);
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
      do_zeros(tot,num_zeros1,q,infile);
      //printf("%lu %lu* %lu",q,index,num_zeros+num_zeros1);
      mpfi_add(tot1,tot1,tot); // add the total twice, one for each char
      mpfi_add(tot1,tot1,tot);
      //mpfi_print_str(" sum ",tot);
    }
  printf("For q= %lu we have sum/phi_2(q) in ",q);
  mpfi_div_ui(tot,tot1,num_prims);
  mpfi_print_str("",tot); 

  return(0);
}

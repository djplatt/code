/*

File: sum_rho.c

Created: Feb 2014

Version: Original

Last Modified: 

Dialect: C

Requires: GMP
          MPFR
          MPFI

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

sum 1/rho+1/(1-rho) and 1/rho

Used for Kadiri data

V 1.0 Initial implementation


By: DJ Platt
    Heilbronn Institute
    Bristol University

Copyright 2014

 */

#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "time.h"
#include "../includes/mpfi_c.h"
#include "../windowed/win_zeta.h"

inline void next_rho(mpfi_t del_t, FILE *infile)
{
  in_bytes(del_t,infile);
}

int main(int argc, char **argv)
{ 
  uint64_t lim;
  if(argc==3)
    lim=atol(argv[2]);
  else if(argc!=2)
    {
      printf("Usage:- %s <zeros file> [t_lim]\n",argv[0]);
      exit(0);
    }

  FILE *infile;
  if(!(infile=fopen(argv[1],"rb")))
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }
  printf("Processing zeros file %s\n",argv[1]);
  mpfi_c_setup(300);
  mpfi_t t;mpfi_init(t);
  mpfi_t del_t;mpfi_init(del_t);
  mpfi_t res;mpfi_init(res);mpfi_set_ui(res,0);
  mpfi_t two;mpfi_init(two);mpfi_set_ui(two,2);
  mpfi_t(tmp);mpfi_init(tmp);
  mpfi_t(res1);mpfi_init(res1);mpfi_set_ui(res1,0);
  mpfi_t (res2);mpfi_init(res2);mpfi_set_ui(res2,0);  

  init_in_bytes();
  //mpfi_print_str("pm1=",pm1);
  unsigned long int num_its;
  fread(&num_its,sizeof(long int),1,infile);
  //printf("Doing %ld iterations.\n",num_its);
  int calced=(1==0);
  double st[2],t0;
  unsigned long int zs[2],it,z;
  int all_done=0;
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile);
      fread(zs,sizeof(long int),1,infile);
      if(st[0]==0.0)
	continue;
      if(!calced)
	{
	  calced=(1==1);
	  t0=st[0];
	}
      fread(zs+1,sizeof(long int),1,infile);
      mpfi_set_d(t,st[0]);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  next_rho(del_t,infile);
          if(mpfi_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  mpfi_add(t,t,del_t); // exact, used to avoid accumulation of +/- 2^-(OP_ACC+1)
	  if((argc==3)&&(mpfi_cmp_ui(t,lim)>0))
	    {printf("Limit reached.\n");all_done=1;break;}
	  mpfi_add(tmp,t,pm1); // not exact
	  mpfi_print_str("gamma   =",tmp);

	  // compute sum 1/gamma
	  mpfi_div(del_t,two,tmp); // 2/gamma
	  mpfi_add(res2,res2,del_t);
	  
	  mpfi_sqr(del_t,tmp);
	  mpfi_add_d(del_t,del_t,0.25);
	  mpfi_div(del_t,two,del_t); // 2/(0.25+t^2)
	  mpfi_add(res,res,del_t);
	  
	}
      if(all_done) break;
    }
  if(argc!=3)
    lim=st[1];
  printf("sum |Im rho| in [%lu,%lu] 1/(rho(1-rho))= ",(uint64_t)t0,lim);  
  mpfi_print_str("",res);
  mpfi_print_str("            1/gamma= ",res2);
  /*
  FILE *outfile=fopen("Ustar_1.dat","wb");
  uint64_t q=1;
  fwrite(&q,sizeof(uint64_t),1,outfile);
  mpfr_t fr;
  mpfr_init(fr);
  mpfi_get_left(fr,res);
  double d=mpfr_get_d(fr,GMP_RNDD);
  fwrite(&d,sizeof(double),1,outfile);
  mpfi_get_right(fr,res);
  d=mpfr_get_d(fr,GMP_RNDU);
  fwrite(&d,sizeof(double),1,outfile);
  fclose(outfile);
  */
  return(0);
}

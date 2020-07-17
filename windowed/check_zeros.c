#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "time.h"
#include "../includes/mpfi_c.h"
#include "../includes/pi_x.h"
#include "../windowed/win_zeta.h"

//#define DEBUG

mpfi_t t,del_t,sig_gap,sig_gap2;

inline void next_rho(mpfi_t del_t, FILE *infile)
{
  in_bytes(del_t,infile);
}


void check_file (char *ifname, FILE *infile)
{
  long int num_blocks,n,zeros[2],z,num_zeros=0;
  double st[2],sig_gap,sig_gap2,norm_gap;
  unsigned char zero[13];
  printf("Checking file %s\n",ifname);
  fread(&num_blocks,sizeof(long int),1,infile);
#ifdef DEBUG
  printf("num blocks=%ld\n",num_blocks);
#endif
  for(n=0;n<num_blocks;n++)
    {
      fread(st,sizeof(double),2,infile);
#ifdef DEBUG
      printf("start=%12.0f finish=%12.0f\n",st[0],st[1]);
#endif
      fread(zeros,sizeof(long int),1,infile);
      mpfi_set_d(t,st[0]);
      if(st[0]!=0)
	{
	  fread(&zeros[1],sizeof(long int),1,infile);
	  num_zeros+=zeros[1]-zeros[0];
	  for(z=zeros[0];z<zeros[1];z++)
	    {
	      //printf("Reading zero number %ld\n",z);
	      next_rho(del_t,infile);
	      /*
	      if(mpfi_cmp_d(del_t,1e-3)<0)
		mpfi_print_str("delta t =",del_t);
	      */
	      mpfi_add(t,t,del_t);
	      //norm_gap=mpfi_get_d(del_t)/(2.0*M_PI)*log(mpfi_get_d(t)/(2.0*M_PI));
	      //sig_gap+=norm_gap;
	      //sig_gap2+=norm_gap*norm_gap;
	    }
	}
#ifdef DEBUG
      mpfi_print_str("Last zero at ",t);
      printf("End of block at %30.28e\n",st[1]);
#endif
      if(mpfi_cmp_d(t,st[1])>0)
	{
	  mpfi_print_str("t=",t);
	  printf("End of block should be %30.28e\n",st[1]);
	}
    }
  printf("Zeros in file %lu\n",num_zeros);
  //printf("sigma gap=%10.8e\n",sig_gap);
  //printf("sigma gap^2=%10.8e\n",sig_gap2);
  //double mean=sig_gap/num_zeros;
  //printf("mean=    %10.8e\n",mean);
  //printf("variance=%10.8e\n",sig_gap2/num_zeros-mean);
}

int main(int argc, char **argv)
{
  FILE *infile;
  mpfr_set_default_prec(200);
  mpfi_init(del_t);
  mpfi_init(t);
  mpfi_init(sig_gap);
  mpfi_init(sig_gap2);
  mpfi_set_ui(sig_gap,0);
  mpfi_set_ui(sig_gap2,0);
  if(argc!=2)
    {
      printf("Usage:- check_zeros <filename>\n");
      exit(0);
    }
  infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("failed to open first file %s. exiting.\n",argv[1]);
      exit(0);
    }
  check_file(argv[1],infile);
  return(0);
}

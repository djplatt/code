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

mpfi_t t,del_t;

inline void next_rho(mpfi_t del_t, FILE *infile)
{
  in_bytes(del_t,infile);
}


void check_file (char *ifname, FILE *infile)
{
  long int num_blocks,n,zeros[2],z,num_zeros=0;
  double st[2],sig_gap,sig_gap2,norm_gap;
  unsigned char zero[13];
  //printf("Checking file %s\n",ifname);
  fread(&num_blocks,sizeof(long int),1,infile);
  for(n=0;n<num_blocks;n++)
    {
      fread(st,sizeof(double),2,infile);
      fread(zeros,sizeof(long int),2,infile);
      mpfi_set_d(t,st[0]);
      num_zeros+=zeros[1]-zeros[0];
      for(z=zeros[0];z<zeros[1];z++)
	{
	  //printf("Reading zero number %ld\n",z);
	  next_rho(del_t,infile);
	  mpfi_add(t,t,del_t);
	  printf("%20.18e\n",mpfi_get_d(t));
	}
    }
}

int main(int argc, char **argv)
{
  FILE *infile;
  mpfr_set_default_prec(200);
  mpfi_init(del_t);
  mpfi_init(t);
  if(argc!=2)
    {
      printf("Usage:- print_zeros <filename>\n");
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

#include "stdlib.h"
#include "stdio.h"
#include "inttypes.h"
#include "../includes/mpfi_c.h"
//#include "gmp.h"
#define PREC (300)

int main(int argc, char** argv)
{
  if(argc!=4)
    {
      printf("Usage:- %s <q> <hurwitz file> <outfile>. Exiting.\n",argv[0]);
      exit(0);
    }
 
  uint64_t q=atol(argv[1]);
  if((q&3)==2)
    {
      printf("moduli congruent to 2 mod 4 have no primitive characters. Exiting.\n");
      exit(0);
    }

  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Failed to open %s for binary input. Exiting.\n",argv[2]);
      exit(0);
    }
  FILE *outfile=fopen(argv[3],"wb");
  if(!outfile)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[3]);
      exit(0);
    }

  uint64_t N_COLUMNS,N_ROWS,M;
  double t_start,t_end,del_t;
  if(fread(&N_COLUMNS,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading number of columns. Exiting.\n");
      exit(0);
    }
  if(fread(&N_ROWS,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading number of rows. Exiting.\n");
      exit(0);
    }
  if(fread(&M,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading M. Exiting.\n");
      exit(0);
    }
  if(fread(&t_start,sizeof(double),1,infile)!=1)
    {
      printf("Error reading t_start. Exiting.\n");
      exit(0);
    }
  if(fread(&t_end,sizeof(double),1,infile)!=1)
    {
      printf("Error reading t_end. Exiting.\n");
      exit(0);
    }
  if(fread(&del_t,sizeof(double),1,infile)!=1)
    {
      printf("Error reading del_t. Exiting.\n");
      exit(0);
    }
  uint64_t num_s=(t_end-t_start)/del_t;
  mpfi_c_setup(PREC);
  mpfi_c_t *hur_vec= (mpfi_c_t *)malloc(sizeof(mpfi_c_t)*N_COLUMNS*N_ROWS*num_s);
  if(!hur_vec)
    {
      printf("Error allocating memory for hur_vec. Exiting.\n");
      exit(0);
    }
  mpz_t l,r;
  mpz_init(l);mpz_init(r);
  for(uint64_t s=0,ptr=0;s<num_s;s++)
    for(uint64_t row=0;row<N_ROWS;row++)
      for(uint64_t column=0;column<N_COLUMNS;column++)
	{
	  if(!mpz_inp_raw(l,infile))
	    {
	      printf("Error reading from hurwitz file. Exiting.\n");
	      exit(0);
	    }
	  if(!mpz_inp_raw(r,infile))
	    {
	      printf("Error reading from hurwitz file. Exiting.\n");
	      exit(0);
	    }
	  mpfi_c_init(hur_vec[ptr]);
	  mpfi_interv_z(hur_vec[ptr]->re,l,r);
	  if(!mpz_inp_raw(l,infile))
	    {
	      printf("Error reading from hurwitz file. Exiting.\n");
	      exit(0);
	    }
	  int64_t exp=mpz_get_si(l);
	  mpfi_mul_2si(hur_vec[ptr]->re,hur_vec[ptr]->re,exp);

	  if(!mpz_inp_raw(l,infile))
	    {
	      printf("Error reading from hurwitz file. Exiting.\n");
	      exit(0);
	    }
	  if(!mpz_inp_raw(r,infile))
	    {
	      printf("Error reading from hurwitz file. Exiting.\n");
	      exit(0);
	    }
	  mpfi_interv_z(hur_vec[ptr]->im,l,r);
	  if(!mpz_inp_raw(l,infile))
	    {
	      printf("Error reading from hurwitz file. Exiting.\n");
	      exit(0);
	    }
	  exp=mpz_get_si(l);
	  mpfi_mul_2si(hur_vec[ptr]->im,hur_vec[ptr]->im,exp);
	  ptr++;
	}
  return(0);
}
	  
	  

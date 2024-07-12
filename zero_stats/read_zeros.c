// read_zeros.c
// example 'C' code to read zeros from lmfdb files and
// convert them into ARB
//
// assuming ARB is installed somewhere sensible
// gcc -O2 read_zeros.c -o read_zeros -larb
// should build it.
//
#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"

// use ARB to manage rounding errors
#include "flint/arb.h"

// the zeros are stored with 101 bits of precision
#define OP_ACC (101)


// the deltas between the centres of intervals for two
// zeros are stored as 13 bytes and are exact
// they are structured as unsigned ints, a:8,b:4,c:1
// the required delta is (c*2^96+b*2^64+a)/2^101
void in_bytes(arb_ptr t, FILE *infile, int64_t prec)
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
  arb_mul_2exp_si(t,t,32); // *2^32
  arb_add_ui(t,t,b,prec);
  arb_mul_2exp_si(t,t,64); // *2^64
  arb_add_ui(t,t,a,prec);
  arb_mul_2exp_si(t,t,-OP_ACC); // *2^(-101)
}

int main(int argc, char ** argv)
{
  printf("Command line was:- ");
  for(unsigned long int i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");

  if(argc!=4)
    {
      printf("Usage:- %s <zeros file list> <N> <prec>\n",argv[0]);
      exit(0);
    }
  
  FILE *lfile=fopen(argv[1],"rb");
  if(lfile == NULL)
    {
      perror("Exiting: ");
      exit(0);
    }
  
  uint64_t N=atol(argv[2]); // how many zeros to read
  arb_t gamma;
  arb_init(gamma);
  int64_t prec=atol(argv[3]); // prec to use within ARB. Sensible to use >> 101

  uint64_t zz=0; // counter for zeros
  char fname[1024];
  arb_t del_t,t; // del_t is gap between zeros. t is centre of zero.
  arb_init(t);arb_init(del_t);
  arb_t z_err;
  arb_init(z_err);
  arb_set_ui(z_err,1);
  arb_mul_2exp_si(z_err,z_err,-OP_ACC-1); // zeros of to within +/- z_err
  arb_t res,tmp,tmp1;
  arb_init(res);arb_init(tmp);arb_init(tmp1);  
  while((zz<N)&&(fscanf(lfile,"%s\n",fname)==1))
    {
      //printf("Processing file %s\n",fname);
      FILE *zfile=fopen(fname,"rb");
      if(zfile == NULL)
	{
	  perror("Exiting: ");
	  exit(0);
	}
      long int num_its;
      // how many blocks does the file contain?
      if(fread(&num_its,sizeof(long int),1,zfile)!=1)
	{
	  printf("Failed to read number of blocks within file %s. Exiting.\n",fname);
	  exit(0);
	}
      double st[2];
      long int zs[2];
      unsigned long int ptr=0;
      
      for(long int it=0;it<num_its;it++)
	{
	  // read start and end of block height
	  if(fread(st,sizeof(double),2,zfile)!=2)
	    {
	      printf("Failed to read start and end of block (doubles). Exiting.\n");
	      exit(0);
	    }
	  // read the number of the first zero in the block
	  if(fread(&zs[0],sizeof(long int),1,zfile)!=1)
	    {
	      printf("Failed to read starting zero count for this block. Exiting.\n");
	      exit(0);
	    }
	  if(st[0]==0.0) // empty block
	    continue;
	  // block contains some zeros, so set start height
	  arb_set_d(t,st[0]);
	  if(zz!=zs[0])
	    {
	      printf("Zero count out of sync. Exiting.\n");
	      exit(0);
	    }
	  // read number of last zero
	  if(fread(&zs[1],sizeof(long int),1,zfile)!=1)
	    {
	      printf("Failed to read ending zero count for this block. Exiting.\n");
	      exit(0);
	    }
	  for(long int z=zs[0]+1;z<=zs[1];z++)
	    {
	      zz++;
	      if(zz>N) // we've read the number of zeros requested
		break;
	      in_bytes(del_t,zfile,prec); // read delta to centre for next zero
	      if(arb_contains_zero(del_t))
		{
		  printf("Two zeros 0 apart. Exiting.\n");
		  exit(0);
		}
	      arb_add(t,t,del_t,prec); // add delta to height
	      arb_set(gamma,t);
	      arb_add_error(gamma,z_err); // add the error for gamma
	    }
	  if(zz>N)
	    break;
	}
      fclose(zfile);
    }
  
  if(zz<N)
    {
      printf("Ran out of zeros after %lu\n",zz);
      exit(0);
    }

  printf("%lu'th is at ",N);
  arb_printd(gamma,40);
  printf("\n");

  return 0;
}


#include "stdio.h"
#include "stdlib.h"
#include "acb.h"
#include "inttypes.h"

//
// Code to demonstrate iterating over the zeros in the LMFDB
// Needs ARB from http://arblib.org/
//

// build it doing something along the lines of
//
// gcc -O2 -o arb_do_zero arb_do_zero.c -larb
//

// the zeros are stored to this accuracy
#define OP_ACC (101)
// in other words, the zeros are accurate to +/- 2^{-102}

// we use a stack to sum things
// this avoids some loss of precision
uint64_t sp,s_count;
#define STACK_LEN (100)
arb_t stack[STACK_LEN];

// call this once at start to set up the stack
void init_stack()
{
  for(uint64_t i=0;i<STACK_LEN;i++)
    arb_init(stack[i]);
  sp=0;
  s_count=0;
}

// call this to add an item into the stack
void add_one(arb_t x, int64_t prec)
{
  if(sp==STACK_LEN)
    {
      printf("Stack overflow. Exiting.\n");
      exit(0);
    }
  arb_set(stack[sp],x);
  s_count++;
  sp++;
  if(!(s_count&1))
    {
      uint64_t s=s_count;
      while(!(s&1))
	{
	  arb_add(stack[sp-2],stack[sp-2],stack[sp-1],prec);
	  sp--;
	  s>>=1;
	}
    }
}

// call this after adding the last item to get the result
void stack_sum(arb_t res, int64_t prec)
{
  arb_set(res,stack[0]);
  for(uint64_t p=1;p<sp;p++)
    arb_add(res,res,stack[p],prec);
}

// read a 13 byte number from file
// structured 8,4,1
// read as if its exact
// represents gap between zeros
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
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (b). Exiting.\n");
      exit(0);
    }
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (c). Exiting.\n");
      exit(0);
    }

  // if I wanted a double, then return (((c*2^32)+b)*2^64+a)/2^101
  // this is the gap from the last zero (or start of file) to this one

  arb_set_ui(t,c);
  arb_mul_2exp_si(t,t,32);
  arb_add_ui(t,t,b,prec);
  arb_mul_2exp_si(t,t,64);
  arb_add_ui(t,t,a,prec);
  arb_mul_2exp_si(t,t,-OP_ACC);
}

// we are going to sum 1/gamma
// put whatever processing you want in here
// it will get called with every gamma in increasing order
void process_zero(arb_t gamma, int64_t prec)
{
  static int init=(1==0);
  static arb_t tmp;
  if(!init)
    {
      init=(0==0);
      arb_init(tmp);
    }
  arb_inv(tmp,gamma,prec); // tmp <- 1/gamma
  add_one(tmp,prec); // put it on the stack
}

int main(int argc, char ** argv)
{
  printf("Command line was:- ");
  for(unsigned long int i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");

  if(argc!=4)
    {
      printf("Usage:- %s <zeros file list> <No of zeros to process, -1 for all> <prec>\n",argv[0]);
      exit(0);
    }
  FILE *lfile=fopen(argv[1],"rb"); // should contain list of zeros files
  if(!lfile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }
  int64_t N=atol(argv[2]); // how many zeros to process, -1 for all
  arb_t T;
  arb_init(T);
  int64_t prec=atol(argv[3]); // working precision. 200 bits is plenty
  uint64_t zz=0;
  char fname[1024];
  arb_t del_t,t;
  arb_init(t);arb_init(del_t);
  arb_t z_err;
  arb_set_ui(z_err,1);
  arb_mul_2exp_si(z_err,z_err,-OP_ACC-1); // this is the error in each zero
  arb_t res,tmp,tmp1;
  arb_init(res);arb_init(tmp);arb_init(tmp1);  
  init_stack();
  // iterate over all files listed in fname
  while(((N<0)||(zz<N))&&(fscanf(lfile,"%s\n",fname)==1))
    {
      FILE *zfile=fopen(fname,"rb");
      if(zfile == NULL)
	{
	  perror("Exiting: ");
	  exit(0);
	}
      uint64_t num_blocks;
      fread(&num_blocks,sizeof(uint64_t),1,zfile); // how many blocks in file
      double st[2];
      uint64_t zs[2];
      // iterate over every block in file
      for(uint64_t block=0;block<num_blocks;block++)
	{
	  fread(st,sizeof(double),2,zfile); // start/end of block as doubles
	  fread(&zs[0],sizeof(uint64_t),1,zfile); // first zero number
	  if(st[0]==0.0) // empty block
	    continue;
	  arb_set_d(t,st[0]);
	  fread(&zs[1],sizeof(uint64_t),1,zfile); // last zero number
	  // iterate over every zero in block
	  for(uint64_t z=zs[0]+1;z<=zs[1];z++)
	    {
	      zz++;
	      if(zz!=z)
		{
		  printf("zz=%lu z=%lu\n",zz,z);
		  exit(0);
		}
	      if((N>=0)&&(zz>N))
		break;
	      in_bytes(del_t,zfile,prec); // gap to the next zero
	      if(arb_contains_zero(del_t))
		{
		  printf("Two zeros 0 apart. Exiting.\n");
		  exit(0);
		}
	      arb_add(t,t,del_t,prec); // this is exact
	      arb_set(tmp,t);
	      arb_add_error(tmp,z_err); // this has +/- on it
	      process_zero(tmp,prec);
	    }
	  if((N>=0)&&(zz>N)) // done all zeros requested
	    break;
	}
      fclose(zfile);
    }

  if((N>=0)&&(zz<N))
    {
      printf("Ran out of zeros after %lu\n",zz);
      exit(0);
    }

  printf("Final zero read was at ");
  arb_printd(tmp,40);
  printf("\n");

  stack_sum(res,prec);
  printf("sum 1/gamma =");arb_printd(res,40);printf("\n");

  return 0;

}


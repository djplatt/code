#include "stdio.h"
#include "stdlib.h"
#include "acb.h"
//#include "../trace/quad.h"

// the zeros are stored with 101 bits of precision
#define OP_ACC (101)

#include "inttypes.h"

// we use a simple stack to sum
uint64_t sp,s_count;
#define STACK_LEN (100)
arb_t *stack;

void init_stack()
{
  stack=(arb_t *)malloc(sizeof(arb_t)*STACK_LEN);
  if(!stack)
    {
      printf("Failed to find memory for stack. Exiting.\n");
      exit(0);
    }
  for(uint64_t i=0;i<STACK_LEN;i++)
    arb_init(stack[i]);
  sp=0;
  s_count=0;
}

// add a new item to the stack
// once we have two items, we add them
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
	  //printf("s_count=%lu adding %lu and %lu\n",s_count,sp-2,sp-1);
	  arb_add(stack[sp-2],stack[sp-2],stack[sp-1],prec);
	  sp--;
	  s>>=1;
	}
    }
}

// add up what's now on the stack.
void stack_sum(arb_t res, int64_t prec)
{
  arb_set(res,stack[0]);
  for(uint64_t p=1;p<sp;p++)
    arb_add(res,res,stack[p],prec);
}



// the deltas between the centres of intervals for two
// zeros are stored in 13 bytes and are exact
// structured as unsigned ints, a:8,b:4,c:1
// the required delta is (c*2^96+b*2^64+c)/2^101
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
  arb_mul_2exp_si(t,t,32);
  arb_add_ui(t,t,b,prec);
  arb_mul_2exp_si(t,t,64);
  arb_add_ui(t,t,a,prec);
  arb_mul_2exp_si(t,t,-OP_ACC);
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
  if(!lfile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }
  uint64_t N=atol(argv[2]);
  arb_t gamma;
  arb_init(gamma);
  int64_t prec=atol(argv[3]);
  uint64_t zz=0;
  char fname[1024];
  arb_t del_t,t;
  arb_init(t);arb_init(del_t);
  arb_t z_err;
  arb_init(z_err);
  arb_set_ui(z_err,1);
  arb_mul_2exp_si(z_err,z_err,-OP_ACC-1);
  arb_t res,tmp,tmp1;
  arb_init(res);arb_init(tmp);arb_init(tmp1);  
  init_stack();
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
      fread(&num_its,sizeof(long int),1,zfile);
      double st[2];
      long int zs[2];
      unsigned long int ptr=0;
      
      for(long int it=0;it<num_its;it++)
	{
	  // read start and end of block height
	  fread(st,sizeof(double),2,zfile);
	  // read the number of the first zero in the block
	  fread(&zs[0],sizeof(long int),1,zfile);
	  if(st[0]==0.0) // empty block
	    continue;
	  // block contains some zeros, so set start height
	  arb_set_d(t,st[0]);
	  // read number of last zero
	  fread(&zs[1],sizeof(long int),1,zfile);
	  for(long int z=zs[0]+1;z<=zs[1];z++)
	    {
	      zz++;
	      if(zz!=z) // these should always be in sync
		{
		  printf("zz=%lu z=%lu\n",zz,z);
		  exit(0);
		}
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
	      // now do something useful with gamma
	      //arb_sqr(tmp1,gamma,prec);
	      //arb_inv(tmp,tmp1,prec); // tmp=1/gamma^2
	      //add_one(tmp,prec);
	    }
	  if(zz>N)
	    break;
	}
      fclose(zfile);
    }

  //printf("zz=%lu\n",zz);
  
  if(zz<N)
    {
      printf("Ran out of zeros after %lu\n",zz);
      exit(0);
    }

  printf("%lu'th is at ",N);
  arb_printd(gamma,40);
  printf("\n");

  //stack_sum(res,prec);

  //printf("sum 1/gamma^2 =");arb_printd(res,40);printf("\n");
  return 0;
}


#include "acb.h"
#include "stdio.h"
#include "inttypes.h"

#define OP_ACC ((int64_t) 101)

// read our 13 byte structure representing a zero gap
// into an int_double
void in_bytes(arb_ptr res, FILE *infile, int64_t prec)
{
  uint64_t a;uint32_t b;uint8_t c;
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
  arb_set_ui(res,c);
  arb_mul_2exp_si(res,res,16);
  arb_add_ui(res,res,b,prec);
  arb_mul_2exp_si(res,res,64);
  arb_add_ui(res,res,a,prec);
  arb_mul_2exp_si(res,res,-OP_ACC);
  return;
}

int main(int argc, char ** argv)
{
  printf("Command line: ");
  for(int64_t i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=3) exit(0);
  FILE*zfile=fopen(argv[1],"rb");
  if(!zfile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }

  int64_t prec=atol(argv[2]);

  long int num_its;
  fread(&num_its,sizeof(long int),1,zfile);
  double st[2];
  long int zs[2];
  int64_t ptr=0;
  arb_t del_t,t,rho;
  arb_init(del_t);arb_init(t);arb_init(rho);arb_init(t_err);
  arb_set_ui(t_err,1);
  arb_mul_2exp_si(t_err,t_err,-OP_ACC-1);
  for(uint64_t it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,zfile);
      fread(&zs[0],sizeof(long int),1,zfile);
      if(st[0]==0.0)
	continue;
      arb_set_d(t,st[0]);
      fread(&zs[1],sizeof(long int),1,zfile);
      for(uint64_t z=zs[0]+1;z<=zs[1];z++)
	{
	  in_bytes(del_t,zfile,prec);
	  if(arb_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  arb_add(t,t,del_t,prec);
	  arb_set(rho,t);
	  arb_add_err(rho,t_err);
	  arb_printd(rho,30);printf("\n");
	}
    }
  return 0;
}

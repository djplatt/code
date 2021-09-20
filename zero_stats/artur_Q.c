#include "stdio.h"
#include "stdlib.h"
#include "acb.h"
//#include "../trace/quad.h"

#define OP_ACC (101)
#include "inttypes.h"

uint64_t sp[2],s_count[2];
#define STACK_LEN (100)
arb_t stack[2][STACK_LEN];

void init_stack(int n)
{
  for(uint64_t i=0;i<STACK_LEN;i++)
    arb_init(stack[n][i]);
  sp[n]=0;
  s_count[n]=0;
}

void add_one(arb_t x, int n, int64_t prec)
{
  if(sp[n]==STACK_LEN)
    {
      printf("Stack overflow. Exiting.\n");
      exit(0);
    }
  arb_set(stack[n][sp[n]],x);
  s_count[n]++;
  sp[n]++;
  if(!(s_count[n]&1))
    {
      uint64_t s=s_count[n];
      while(!(s&1))
	{
	  //printf("s_count=%lu adding %lu and %lu\n",s_count,sp-2,sp-1);
	  arb_add(stack[n][sp[n]-2],stack[n][sp[n]-2],stack[n][sp[n]-1],prec);
	  sp[n]--;
	  s>>=1;
	}
    }
}

void stack_sum(arb_t res, int n, int64_t prec)
{
  arb_set(res,stack[n][0]);
  for(uint64_t p=1;p<sp[n];p++)
    arb_add(res,res,stack[n][p],prec);
}

// read a 13 byte number from file
// structured 8,4,1
// read as if its exact
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

void do_gamma(arb_t res, arb_t gam, arb_t s, int64_t prec)
{
  static int init=(1==0);
  static arb_t tmp1,tmp2,tmp3,tmp4,qtr;
  if(!init)
    {
      init=(1==1);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(tmp4);
      arb_init(qtr);
      arb_set_d(qtr,0.25);
    }
  arb_sqr(tmp1,gam,prec);
  arb_add(tmp2,tmp1,qtr,prec);
  arb_mul_2exp_si(tmp1,s,-1); // s/2
  arb_neg(tmp1,tmp1); // -s/2
  arb_pow(tmp3,tmp2,tmp1,prec); //(1/4+gam^2)^-s/2
  arb_mul_2exp_si(tmp1,gam,1); // 2gam
  arb_atan(tmp2,tmp1,prec); // atan(2gam)
  arb_mul(tmp4,tmp2,s,prec); // s atan (2 gam)
  arb_sin(tmp1,tmp4,prec); // sin
  arb_mul(res,tmp1,tmp3,prec);
  //add_one(tmp2,0,prec);
}

int main(int argc, char ** argv)
{
  printf("Command line was:- ");
  for(unsigned long int i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");

  if(argc!=6)
    {
      printf("Usage:- %s <zeros file list> <N> <sn> <sd> <prec>\n",argv[0]);
      exit(0);
    }
  FILE *lfile=fopen(argv[1],"rb");
  if(!lfile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }
  uint64_t N=atol(argv[2]);
  arb_t T;
  arb_init(T);
  int64_t prec=atol(argv[5]);
  arb_t s;
  arb_init(s);
  arb_set_ui(s,atol(argv[3]));
  arb_div_ui(s,s,atol(argv[4]),prec);
  printf("s set to ");arb_printd(s,10);printf("\n");
  uint64_t zz=0;
  char fname[1024];
  arb_t del_t,t;
  arb_init(t);arb_init(del_t);
  arb_t z_err;
  arb_set_ui(z_err,1);
  arb_mul_2exp_si(z_err,z_err,-OP_ACC-1);
  arb_t res1,res2,tmp,tmp1;
  arb_init(res1);arb_init(res2);arb_init(tmp);arb_init(tmp1);  
  init_stack(0);init_stack(1);
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
      fread(&num_its,sizeof(long int),1,zfile);
      double st[2];
      long int zs[2];
      unsigned long int ptr=0;
      
      for(long int it=0;it<num_its;it++)
	{
	  //if((it%100)==0)
	  //printf("Starting block %lu/%lu\n",it+1,num_its);
	  fread(st,sizeof(double),2,zfile);
	  fread(&zs[0],sizeof(long int),1,zfile);
	  if(st[0]==0.0)
	    continue;
	  arb_set_d(t,st[0]);
	  fread(&zs[1],sizeof(long int),1,zfile);
	  //printf("Processing zero %ld to %ld=%ld in total.\n",zs[0]+1,zs[1],zs[1]-zs[0]);
	  for(long int z=zs[0]+1;z<=zs[1];z++)
	    {
	      zz++;
	      if(zz!=z)
		{
		  printf("zz=%lu z=%lu\n",zz,z);
		  exit(0);
		}
	      if(zz>N)
		break;
	      in_bytes(del_t,zfile,prec);
	      if(arb_contains_zero(del_t))
		{
		  printf("Two zeros 0 apart. Exiting.\n");
		  exit(0);
		}
	      arb_add(t,t,del_t,prec);
	      arb_set(tmp,t);
	      arb_add_error(tmp,z_err);
	      do_gamma(res2,tmp,s,prec);
	      arb_add(res1,res1,res2,prec);
	      arb_sub_ui(res2,res1,1,prec);
	      if(arb_is_positive(res2))
		{
		  printf("Partial sum exceed 1 at zero %lu ",z);
		  arb_printd(res1,20);
		  printf("\n");
		  exit(0);
		}
	      if(!arb_is_negative(res2))
		{
		  printf("Partial sum straddled 1. Exiting.\n");
		  return 0;
		}

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

  arb_set(T,t);
  printf("%lu'th is at ",N);
  arb_printd(T,40);
  printf("\n");
  printf("Sum so far = ");
  arb_printd(res1,20);
  printf("\n");

  //stack_sum(res1,0,prec);
  //arb_mul_2exp_si(res1,res1,1);
  //stack_sum(res2,1,prec);
  //arb_neg(res2,res2);

  //printf("sum  =");arb_printd(res1,40);printf("\n");
  //printf("sum Im 1/(1/2+igam)^s =");arb_printd(res2,40);printf("\n");
  return 0;
}


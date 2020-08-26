#include "stdio.h"
#include "stdlib.h"
#include "acb_hypgeom.h"
//#include "../trace/quad.h"

#define OP_ACC (101)
#include "inttypes.h"

uint64_t sp,s_count;
#define STACK_LEN (100)
arb_t stack[STACK_LEN];

void init_stack()
{
  for(uint64_t i=0;i<STACK_LEN;i++)
    arb_init(stack[i]);
  sp=0;
  s_count=0;
}

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

void stack_sum(arb_t res, int64_t prec)
{
  arb_set(res,stack[0]);
  for(uint64_t p=1;p<sp;p++)
    arb_add(res,res,stack[p],prec);
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



arb_t two_pi;
void E2(arb_t res, arb_t T, int64_t prec)
{
  static int init=(1==0);
  static arb_t tmp,tmp1;
  if(!init)
    {
      init=(1==1);
      arb_init(tmp);
      arb_init(tmp1);
    }
  arb_set_ui(tmp,236);
  arb_div_ui(tmp1,tmp,1000,prec);
  arb_log(tmp,T,prec);
  arb_mul(res,tmp,tmp1,prec);
  arb_set_ui(tmp,8527);
  arb_div_ui(tmp1,tmp,1000,prec);
  arb_add(res,res,tmp1,prec);
  arb_div(res,res,T,prec);
  arb_div(res,res,T,prec);
  arb_div(res,res,T,prec);
}  

void Q(arb_t res, arb_t T, uint64_t N,int64_t prec)
{
  static int init=(1==0);
  static arb_t tmp,tmp1;
  if(!init)
    {
      init=(1==1);
      arb_init(tmp);
      arb_init(tmp1);
    }
  arb_set_ui(res,N);
  arb_div(tmp,T,two_pi,prec);
  arb_log(tmp1,tmp,prec);
  arb_sub_ui(tmp1,tmp1,1,prec);
  arb_mul(tmp,tmp,tmp1,prec);
  arb_set_d(tmp1,7.0/8.0);
  arb_add(tmp1,tmp1,tmp,prec);
  arb_sub(res,res,tmp1,prec);
}

void phi(arb_t res, arb_t t, int64_t prec)
{
  static int init=(1==0);
  static arb_t tmp,tmp1;
  if(!init)
    {
      init=(1==1);
      arb_init(tmp);
      arb_init(tmp1);
    }
  arb_sqr(tmp,t,prec);
  arb_inv(res,tmp,prec);
}

void int_fix(arb_t res, arb_t T, int64_t prec)
{
  static int init=(1==0);
  static arb_t tmp,tmp1;
  if(!init)
    {
      init=(1==1);
      arb_init(tmp);
      arb_init(tmp1);
    }
  arb_div(tmp,T,two_pi,prec);
  arb_log(tmp1,tmp,prec);
  arb_add_ui(tmp,tmp1,1,prec);
  arb_div(tmp1,tmp,two_pi,prec);
  arb_div(res,tmp1,T,prec);
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
  arb_t T;
  arb_init(T);
  int64_t prec=atol(argv[3]);
  uint64_t zz=0;
  char fname[1024];
  arb_t del_t,t;
  arb_init(t);arb_init(del_t);
  arb_t z_err;
  arb_set_ui(z_err,1);
  arb_mul_2exp_si(z_err,z_err,-OP_ACC-1);
  arb_t res,tmp,tmp1;
  arb_init(res);arb_init(tmp);arb_init(tmp1);
  arb_init(two_pi);
  arb_const_pi(two_pi,prec);
  arb_mul_2exp_si(two_pi,two_pi,1);
  init_stack();
  FILE *zfile;
  while((zz<N)&&(fscanf(lfile,"%s\n",fname)==1))
    {
      //printf("Processing file %s\n",fname);
      zfile=fopen(fname,"rb");
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
	      phi(tmp1,tmp,prec);
	      add_one(tmp1,prec);
	    }
	  if(zz>N)
	    break;
	}
      if(zz<=N) fclose(zfile);
    }

  //printf("zz=%lu\n",zz);
  
  if(zz<N)
    {
      printf("Ran out of zeros after %lu\n",zz);
      exit(0);
    }

  arb_set(T,t);
  arb_set(tmp,T);
  arb_add_error(tmp,z_err);
  printf("%lu'th zero is at ",N);
  arb_printd(tmp,40);
  printf("\n");
  in_bytes(del_t,zfile,prec); // will fail if last zero in block read last
  fclose(zfile);
  arb_mul_2exp_si(del_t,del_t,-1); // 1/2 way to next zero
  arb_add(T,T,del_t,prec);
  printf("T set to ");arb_printd(T,40);printf("\n");

  stack_sum(res,prec);

  printf("sum 1/gamma^2 =");arb_printd(res,40);printf("\n");
  Q(tmp,T,N,prec);
  printf("Q(T) is ");arb_printd(tmp,40);printf("\n");
  phi(tmp1,T,prec);
  printf("phi(T) is ");arb_printd(tmp1,40);printf("\n");
  arb_mul(tmp,tmp,tmp1,prec);
  printf("phi(T)Q(T) = ");arb_printd(tmp,40);printf("\n");
  arb_sub(res,res,tmp,prec);

  printf("sum 1.gamma^2-phi(T)Q(T) = ");arb_printd(res,40);printf("\n");
  int_fix(tmp,T,prec);
  printf("Integral term = ");arb_printd(tmp,40);printf("\n");
  arb_add(res,res,tmp,prec);
  printf("sum 1.gamma^2-phi(T)Q(T)+1/2pi int phi(t)log(t/2pi) dt = ");arb_printd(res,40);printf("\n");

  E2(tmp,T,prec);
  printf("|E2| = ");arb_printd(tmp,4);printf("\n");
  arb_add_error(res,tmp);
  printf("Final interval is ");arb_printd(res,40);printf("\n");
  return 0;

}


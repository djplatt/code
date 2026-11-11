#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"
#include "flint/acb.h"
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
#define STACK_LEN (100)
#define NUM_STACKS (6)
uint64_t sp[NUM_STACKS],s_count[NUM_STACKS];
arb_t stack[NUM_STACKS][STACK_LEN];

// call this once at start to set up the stack
void init_stacks()
{
  for(uint64_t s=0;s<NUM_STACKS;s++)
    {
      for(uint64_t i=0;i<STACK_LEN;i++)
	arb_init(stack[s][i]);
      sp[s]=0;
      s_count[s]=0;
    }
}

// call this to add an item into the stack
void add_one(arb_t x, uint64_t st, int64_t prec)
{
  if(sp[st]==STACK_LEN)
    {
      printf("Stack overflow. Exiting.\n");
      exit(0);
    }
  arb_set(stack[st][sp[st]],x);
  s_count[st]++;
  sp[st]++;
  if(!(s_count[st]&1))
    {
      uint64_t s=s_count[st];
      while(!(s&1))
	{
	  arb_add(stack[st][sp[st]-2],stack[st][sp[st]-2],stack[st][sp[st]-1],prec);
	  sp[st]--;
	  s>>=1;
	}
    }
}

// call this after adding the last item to get the result
void stack_sum(arb_t res, uint64_t st, int64_t prec)
{
  arb_set(res,stack[st][0]);
  for(uint64_t p=1;p<sp[st];p++)
    arb_add(res,res,stack[st][p],prec);
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
void process_zero(arb_t gamma, acb_t residue, int64_t prec)
{
  static int init=false;
  static arb_t tmp, abs_rho,abs_residue,abs_zeta;
  static arb_t two_k,two_k_zeta,rho_zetap,rho_zetap_zeta;
  static acb_t rho,two_rho,zeta_two_rho;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(abs_rho);
      arb_init(abs_residue);
      arb_init(abs_zeta);
      arb_init(two_k);
      arb_init(two_k_zeta);
      arb_init(rho_zetap);
      arb_init(rho_zetap_zeta);
      acb_init(rho);
      arb_set_d(acb_realref(rho),0.5);
      acb_init(two_rho);
      arb_set_ui(acb_realref(two_rho),1);
      acb_init(zeta_two_rho);
    }
  arb_set(acb_imagref(rho),gamma);
  arb_mul_2exp_si(acb_imagref(two_rho),gamma,1);
  acb_zeta(zeta_two_rho,two_rho,130); // want 100 bits at 10^9
  acb_abs(abs_zeta,zeta_two_rho,prec);
  acb_abs(abs_rho,rho,prec); // |rho|
  acb_abs(abs_residue,residue,prec); // |1/zeta'(rho)|
  arb_div(rho_zetap,abs_residue,abs_rho,prec); // |1/rho zeta'(rho)|
  arb_mul(rho_zetap_zeta,rho_zetap,abs_zeta,prec); // |zeta(2 rho)/(rho zeta'(rho)|
  arb_sqr(two_k,rho_zetap,prec); // |1/...| ^2
  arb_sqr(two_k_zeta,rho_zetap_zeta,prec); // |zeta(2 rho)/...|^2
  add_one(two_k,0,prec);
  add_one(two_k_zeta,3,prec);
  arb_sqr(rho_zetap,two_k,prec); // ^4
  add_one(rho_zetap,1,prec);
  arb_sqr(rho_zetap_zeta,two_k_zeta,prec); // ^4
  add_one(rho_zetap_zeta,4,prec);
  arb_mul(abs_rho,rho_zetap,two_k,prec); // ^6
  add_one(abs_rho,2,prec);
  arb_mul(abs_rho,rho_zetap_zeta,two_k_zeta,prec);
  add_one(abs_rho,5,prec);
}

int main(int argc, char ** argv)
{
  printf("Command line was:- ");
  for(unsigned long int i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");

  if(argc!=4)
    {
      printf("Usage:- %s <residues file list> <No of zeros to process, -1 for all> <prec>\n",argv[0]);
      exit(0);
    }
  FILE *infile=fopen(argv[1],"rb"); // should contain list of zeros files
  if(!infile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }
  int64_t N=atol(argv[2]); // how many zeros to process, -1 for all
  arb_t T;
  arb_init(T);
  int64_t prec=atol(argv[3]); // working precision. 200 bits is plenty
  char fname[1024];
  arb_t res,gamma,tmp1;
  arb_init(res);arb_init(gamma);arb_init(tmp1);  
  init_stacks();
  // iterate over all files listed in fname
  acb_t residue;
  acb_init(residue);
  uint64_t zz=0;
  while(fscanf(infile,"%s\n",fname)==1)
    {
      FILE *rfile=fopen(fname,"rb");
      if(rfile == NULL)
	{
	  perror("Exiting: ");
	  exit(0);
	}
      while(arb_load_file(gamma,rfile)==0)
	{
	  if(arb_load_file(acb_realref(residue),rfile)!=0)
	    {
	      printf("Error reading real part of residue. Exiting.\n");
	      return 0;
	    }
	  if(arb_load_file(acb_imagref(residue),rfile)!=0)
	    {
	      printf("Error reading imaginary part of residue. Exiting.\n");
	      return 0;
	    }
	  process_zero(gamma,residue,prec);
	  zz++;
	  if((N>0)&&(zz>=N))
	    {
	      fclose(rfile);
	      break;
	    }
	}
      if((N>0)&&(zz>=N))
	break;
    }
  fclose(infile);
  printf("Final zero read was at ");
  arb_printd(gamma,40);
  printf("\n");

  for(uint64_t i=0;i<6;i++)
    {
      stack_sum(res,i,prec);
      printf("sum %lu = ",i);arb_printd(res,40);printf("\n%lu ",i);
      arb_dump_file(stdout,res);
      printf("\n");
    }

  
  return 0;

}


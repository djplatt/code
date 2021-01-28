#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"
#include "acb.h"
#include "inttypes.h"

//
// look for extremes of S(t)
//

// build it doing something along the lines of
//
// gcc -O2 -o S S.c -larb
//

// the zeros are stored to this accuracy
#define OP_ACC (101)
// in other words, the zeros are accurate to +/- 2^{-102}


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

void arb_S(arb_t res, arb_t nt, arb_t gamma, int64_t prec)
{

  static bool init=false;
  static acb_t z,z1;
  static arb_t pi,lnpi,tmp1,tmp2;
  if(!init)
    {
      init=true;
      acb_init(z);
      acb_init(z1);
      arb_set_d(acb_realref(z),0.25);
      arb_init(pi);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_const_pi(pi,prec);
      arb_init(lnpi);
      arb_log(lnpi,pi,prec);
    }
  //printf("In arb_S with gamma = ");arb_printd(gamma,20); printf(" N(t) = ");arb_printd(nt,20);printf("\n");
  arb_set(acb_imagref(z),gamma);
  arb_mul_2exp_si(acb_imagref(z),acb_imagref(z),-1); // 1/4+it/2
  acb_lgamma(z1,z,prec);
  arb_mul(tmp1,acb_imagref(z),lnpi,prec);
  arb_sub(tmp2,acb_imagref(z1),tmp1,prec);
  arb_div(tmp1,tmp2,pi,prec);
  arb_sub(tmp2,nt,tmp1,prec);
  arb_sub_ui(res,tmp2,1,prec);
}
  
// it will get called with every gamma in increasing order
bool process_zero(arb_t gamma, uint64_t z, arb_t max, int64_t prec)
{
  static bool init=false;
  static arb_t tmp,tmp1,tmp2,tmp3,nt;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(nt);
    }
  arb_set_ui(nt,z);
  //printf("In process_zero with gamma = ");arb_printd(gamma,20); printf(" N(t) = ");arb_printd(nt,20);printf("\n");
  arb_S(tmp,nt,gamma,prec);
  //printf("S(t) = ");arb_printd(tmp1,20);printf("\n");
  arb_abs(tmp1,tmp);
  arb_sub_ui(tmp2,tmp,1,prec);
  arb_abs(tmp3,tmp2);
  arb_max(tmp2,tmp1,tmp3,prec);
  arb_sub(tmp1,max,tmp2,prec);
  if(arb_is_negative(tmp1))
    {
      //printf("New max ");arb_printd(tmp2,20);printf(" at ");arb_printd(gamma,20);printf("\n");
      //printf("S(t) = ");arb_printd(tmp,20);printf("\n");
      arb_set(max,tmp2);
      return true;
    }
  if(arb_is_positive(tmp1))
    return false;
  printf("Error at gamma = ");arb_printd(gamma,20);printf(" can't distinguish S(t) and max.\n");
  exit(0);
  return false;
}

int main(int argc, char ** argv)
{
  printf("Command line was:- ");
  for(unsigned long int i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");

  if(argc!=3)
    {
      printf("Usage:- %s <zero file> <prec>\n",argv[0]);
      exit(0);
    }
  FILE *zfile=fopen(argv[1],"rb"); // should contain list of zeros
  if(!zfile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }
  int64_t prec=atol(argv[2]); // working precision. 200 bits is plenty
  uint64_t zz=0;
  arb_t del_t,t;
  arb_init(t);arb_init(del_t);
  arb_t z_err;
  arb_init(z_err);
  arb_set_ui(z_err,1);
  arb_mul_2exp_si(z_err,z_err,-OP_ACC-1); // this is the error in each zero
  arb_t res,tmp,tmp1;
  arb_init(res);arb_init(tmp);arb_init(tmp1);  
  arb_t max,max_t;
  arb_init(max);
  arb_set_ui(max,1); // S(0)=-1
  arb_init(max_t);
  uint64_t max_z;
  
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
	  in_bytes(del_t,zfile,prec); // gap to the next zero
	  if(arb_contains_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  arb_add(t,t,del_t,prec); // this is exact
	  arb_set(tmp,t);
	  arb_add_error(tmp,z_err); // this has +/- on it
	  if(process_zero(tmp,z,max,prec))
	    {
	      arb_set(max_t,t);
	      max_z=z;
	    }
	}
    }

  printf("Final zero read was at ");
  arb_printd(tmp,40);
  printf("\n");

  printf("Max S(t) seen was ");arb_printd(max,20);printf(" at zero %lu near ",max_z);arb_printd(max_t,20);printf("\n");


  return 0;

}


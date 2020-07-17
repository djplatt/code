// test speed for Chris. See 4/10/16 email
#include "stdio.h"
#include "stdlib.h"
#include "acb.h"
#include "../trace/quad.h"

#define OP_ACC (101)
#include "inttypes.h"

arb_t pm1; // +/- 1 in the last place

// initialise some variables
void init_in_bytes()
{
  arb_t tmp;
  arb_init(tmp);
  arb_init(pm1);
  arb_set_ui(tmp,1);
  arb_mul_2exp_si(tmp,tmp,-OP_ACC-1);
  arb_zero(pm1);
  arb_add_error(pm1,tmp);
  arb_clear(tmp);
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


// read the next imaginary part of rho into del_t
inline void next_rho(arb_t del_t, FILE *infile, int64_t prec)
{
  in_bytes(del_t,infile,prec);
}


// Chris Smith (York) test version
inline void do_rho(arb_t res1, arb_t gamma, arb_t omega, arb_t two_alpha, int64_t prec)
{
  // set up static variables
  static int init=(1==0);
  static arb_t g2,gam_om,tmp,tmp1,tmp2,tmp3,tmp4,gsgo,cgo,den,deninc;
  if(!init)
    {
      init=(1==1);
      arb_init(g2);
      arb_init(gam_om);
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(tmp4);
      arb_init(gsgo);
      arb_init(cgo);
      arb_init(den);
      arb_init(deninc);
    }

  //printf("In do_rho with gamma=");arb_printd(gamma,20);printf("\n");
  arb_mul(g2,gamma,gamma,prec); // gamma^2
  arb_mul(tmp,g2,two_alpha,prec);
  arb_exp(tmp1,tmp,prec); // exp(-gam^2/(2 alpha))
  arb_mul_2exp_si(tmp,g2,2);
  arb_add_ui(tmp2,tmp,1,prec);
  arb_mul_2exp_si(tmp2,tmp2,-2); // (1/4+gamma^2)
  arb_mul(tmp,tmp2,tmp2,prec);  
  arb_mul(tmp2,tmp,omega,prec);
  arb_div(tmp,tmp1,tmp2,prec); // exp(-gamma^2/(2alpha))/( omega (1/4+gamma^2)^2)
  //printf("exp()/blah=");arb_printd(tmp,20);printf("\n");

  arb_mul(gam_om,gamma,omega,prec); // gamma*omega
  arb_sin_cos(tmp1,cgo,gam_om,prec);
  arb_mul(gsgo,tmp1,gamma,prec); // gamma*sin(gamma*omega)
  arb_mul_2exp_si(gsgo,gsgo,1); // 2 gamma*sin(gamma*omega)

  arb_mul_2exp_si(tmp1,g2,2); // 4 gamma^2
  arb_sub_ui(tmp2,tmp1,1,prec);
  arb_mul_2exp_si(tmp2,tmp2,-1); // 2gamma^2-1/2
  arb_mul(tmp1,tmp2,cgo,prec); // (2gamma^2-1/2)cos(omega gamma)
  arb_sub(tmp2,gsgo,tmp1,prec); // (1/2-2gamma^2)cos() + 2gamma sin()

  arb_mul(res1,tmp2,tmp,prec);
  arb_neg(res1,res1);
  //printf("Result=");arb_printd(res1,20);printf("\n");
  //exit(0);
}

int main(int argc, char **argv)
{
  int i;
  int64_t prec=200;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=6)
    {
      printf("Fatal error in main: usage skewes <prec> <zeros file> <int omega> <omega div> <alpha>. Exiting.\n");
      exit(0);
    }
  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Fatal error in main: failed to open zeros file for binary input. Exiting.\n");
      exit(0);
    }

  arb_t gamma,omega,res1,tot1,del_t,t,two_alpha;
  arb_init(gamma);arb_init(omega);arb_init(del_t);arb_init(t);
  arb_init(res1);arb_init(tot1);



  arb_set_ui(omega,atol(argv[3]));
  arb_div_ui(omega,omega,atol(argv[4]),prec);
  //arb_out_str(stdout,10,0,omega);printf("\n");
  
  double alpha_d=atof(argv[5]);//66.0e11;//5444646098003630.0;


  printf("omega = ");arb_printd(omega,20);printf("\n");

  arb_init(two_alpha);

  arb_set_d(two_alpha,-alpha_d);
  printf("-alpha = ");arb_printd(two_alpha,20);printf("\n");
  arb_mul_2exp_si(two_alpha,two_alpha,1);
  arb_inv(two_alpha,two_alpha,prec); // -1/(2*alpha)
  
  init_in_bytes();

  arb_set_ui(tot1,0);

  long int num_its,it,z;
  double st[2];
  long int zs[2],n_zeros=0;
  fread(&num_its,sizeof(long int),1,infile);
  //printf("Doing %ld iterations.\n",num_its);
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile); // starting/ending t, exact
      fread(&zs[0],sizeof(long int),1,infile); // starting zero number
      if(st[0]==0.0)
	continue;
      fread(&zs[1],sizeof(long int),1,infile); // ending zero number
      n_zeros+=zs[1]-zs[0];
      arb_set_d(gamma,st[0]);
      arb_set(t,gamma);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  next_rho(del_t,infile,prec); // distance to next gamma
          if(arb_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  arb_add(t,t,del_t,prec); // exact, used to avoid accumulation of +/- 2^-(OP_ACC+1)
	  arb_add(gamma,t,pm1,prec); // not exact

	  // compute Chris's sum at this gamma
	  do_rho(res1,gamma,omega,two_alpha,prec);
	  // add the results to our running totals
	  arb_add(tot1,tot1,res1,prec);
	  
	}
    }
  printf("We processed %ld zeros.\n",n_zeros);
  printf("test on file %s at omega = %f returned ",(double)atol(argv[3])/atol(argv[4]),argv[2]);
  arb_printd(tot1,prec);
  printf("\n");
  return(0);
}

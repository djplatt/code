/*
  Read gammas and residues from list of files.
  sum all into sum1,sum2,sum3
  if gamma < t0/2
     then sum into sum1a,sum2a and sum3a as well
  t0 should be end of largest file
*/

#include "flint/acb.h"
#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"

#include "inttypes.h"

void eta(arb_t res, arb_t t, int64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_neg(tmp,t);
  arb_add_ui(res,tmp,2,prec);
}

void et_gam(arb_t res, arb_t gamma, uint64_t t0_by_2, int64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_div_ui(tmp,gamma,t0_by_2,prec);
  eta(res,tmp,prec);
}
  
int main(int argc, char **argv)
{
  uint64_t i;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=4)
    {
      printf("Fatal error in %s: usage %s <prec> <file list> <t0>. Exiting.\n",argv[0],argv[0]);
      exit(0);
    }
  int64_t prec=atol(argv[1]);
  printf("ARB working precision set to %ld\n",prec);
  FILE *infile1=fopen(argv[2],"r");
  if(!infile1)
    {
      printf("Fatal error in %s: failed to open list of files %s input. Exiting.\n",argv[0],argv[2]);
      exit(0);
    }
  int64_t t0=atol(argv[3]);
  acb_t res,rho;
  arb_t tmp0,tmp,tmp1,tmp2,tmp3,tmp4,gamma,sum1,sum2,sum3,sum1a,sum2a,sum3a;
  arb_t sum3b,sum2b,sum1b;
  uint64_t zero_count=0;
  acb_init(res);acb_init(rho);
  arb_set_d(acb_realref(rho),0.5); // will be 1/2+i gamma
  arb_init(tmp0);arb_init(tmp);arb_init(tmp1);arb_init(tmp2);arb_init(gamma);
  arb_init(tmp3);arb_init(tmp4);
  arb_init(sum1);arb_init(sum2);arb_init(sum3);
  arb_init(sum1a);arb_init(sum2a);arb_init(sum3a);
  arb_init(sum1b);arb_init(sum2b);arb_init(sum3b);

  FILE *infile;
  char fname[200];
  //uint64_t next=100000; // print sums out every so often
  //uint64_t del=100000;
  bool under; // true if gamma < t0/2.0
  int cusp; // cusp <0 means all gamma < t0/2
            // cusp == 0 means cross over
            // cusp > 0 means all gamma > t0/2
  while(fscanf(infile1,"%d %s\n",&cusp,fname)==2)
    {
      under=true;
      if(cusp>0)
	under=false;

      if(cusp<0) printf("summing fname=%s into sum_i and sum_ia\n",fname);
      if(cusp>0) printf("summing fname=%s into sum_i and sum_ib\n",fname);
      if(cusp==0) printf("summing fname=%s into sum_i and sum_ia to start\n",fname);
      
      infile=fopen(fname,"r");
      if(!infile)
	{
	  printf("Fatal error in %s: failed to open zeros file at %s for input. Exiting.\n",argv[0],fname);
	  exit(0);
	}
      
      while((arb_load_file(gamma,infile)==0)&&(arb_load_file(acb_realref(res),infile)==0)&&(arb_load_file(acb_imagref(res),infile)==0))
	{
	  zero_count++;

	  if(cusp==0)
	    if(under) // only need to see gamma > t0/2 once
	      {
		arb_sub_ui(tmp,gamma,t0/2,prec);
		under=arb_is_negative(tmp); // gamma is still < t0/2
		if(!under)
		  {printf("switched summing from a to b at gamma = ");arb_printd(gamma,20);printf("\n");}
	      }
	  arb_set(acb_imagref(rho),gamma);
	  acb_abs(tmp0,rho,prec);
	  arb_inv(tmp,tmp0,prec);
	  acb_abs(tmp1,res,prec);
	  arb_mul(tmp2,tmp,tmp1,prec); // |res|/|rho|
	  arb_add(sum1,sum1,tmp2,prec);
	  if(under) // gamma < t0/2
	    arb_add(sum1a,sum1a,tmp2,prec);
	  else
	    {
	      et_gam(tmp3,gamma,t0/2,prec);
	      arb_mul(tmp4,tmp2,tmp3,prec);
	      arb_add(sum1b,sum1b,tmp4,prec);
	    }
	    
	  arb_mul(tmp1,tmp2,tmp,prec); // |res|/|rho|^2
	  arb_add(sum2,sum2,tmp1,prec);
	  if(under)
	    arb_add(sum2a,sum2a,tmp1,prec);
	  else
	    {
	      arb_mul(tmp4,tmp1,tmp3,prec);
	      arb_add(sum2b,sum2b,tmp4,prec);
	    }
	  arb_mul(tmp2,tmp1,tmp,prec); // |res|/|rho|^3
	  arb_add(sum3,sum3,tmp2,prec);
	  if(under)
	    arb_add(sum3a,sum3a,tmp2,prec);
	  else
	    {
	      arb_mul(tmp4,tmp2,tmp3,prec);
	      arb_add(sum3b,sum3b,tmp4,prec);
	    }
	}
      
      fclose(infile);
      printf("We processed %lu zeros in %s.\n",zero_count,fname);
      
    }

  fclose(infile1);
  //arb_printd(gamma,20);printf("\n");
  printf("sum 1 ");arb_printd(sum1,20);printf("\nsum 2 ");
  arb_printd(sum2,20);printf("\nsum 3 ");
  arb_printd(sum3,20);printf("\n");
  printf("sum 1a ");arb_printd(sum1a,20);printf("\nsum 2a ");
  arb_printd(sum2a,20);printf("\nsum 3a ");
  arb_printd(sum3a,20);printf("\n");
  printf("sum 1b ");arb_printd(sum1b,20);printf("\nsum 2b ");
  arb_printd(sum2b,20);printf("\nsum 3b ");
  arb_printd(sum3b,20);printf("\n");

  acb_clear(res);
  arb_clear(tmp);
  arb_clear(tmp1);
  arb_clear(tmp2);
  arb_clear(sum1);
  arb_clear(sum2);
  arb_clear(sum3);
  arb_clear(gamma);
  return 0;
}



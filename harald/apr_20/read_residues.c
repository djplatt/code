#include "acb.h"
#include "stdio.h"
#include "stdlib.h"

#include "inttypes.h"

int main(int argc, char **argv)
{
  uint64_t i;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=3)
    {
      printf("Fatal error in %s: usage %s <prec> <file list>. Exiting.\n",argv[0],argv[0]);
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
  acb_t res;
  arb_t tmp,tmp1,tmp2,gamma,sum1,sum2,sum3;
  uint64_t zero_count=0;
  acb_init(res);
  arb_init(tmp);arb_init(tmp1);arb_init(tmp2);arb_init(gamma);
  arb_init(sum1);arb_init(sum2);arb_init(sum3);

  FILE *infile;
  char fname[200];
  uint64_t next=100000; // print sums out every so often
  uint64_t del=100000;
  while(fscanf(infile1,"%s\n",fname)==1)
    {
      
      printf("fname=%s next print point = %lu\n",fname,next);
      infile=fopen(fname,"r");
      if(!infile)
	{
	  printf("Fatal error in %s: failed to open zeros file at %s for input. Exiting.\n",argv[0],fname);
	  exit(0);
	}
      
      while((arb_load_file(gamma,infile)==0)&&(arb_load_file(acb_realref(res),infile)==0)&&(arb_load_file(acb_imagref(res),infile)==0))
	{
	  
	  arb_sub_ui(tmp,gamma,next,prec);
	  if(!arb_is_negative(tmp))
	    {
	      printf("Values at %lu ",next);
	      arb_printd(gamma,20);printf(" ");
	      arb_printd(sum1,20);printf(" ");
	      arb_printd(sum2,20);printf(" ");
	      arb_printd(sum3,20);printf("\n");
	      next+=del;
	    }
	  
	  zero_count++;
	  
	  arb_inv(tmp,gamma,prec);
	  acb_abs(tmp1,res,prec);
	  arb_mul(tmp2,tmp,tmp1,prec); // res/gamma
	  arb_add(sum1,sum1,tmp2,prec);
	  arb_mul(tmp1,tmp2,tmp,prec); // res/gamma^2
	  arb_add(sum2,sum2,tmp1,prec);
	  arb_mul(tmp2,tmp1,tmp,prec); // res/gamma^3
	  arb_add(sum3,sum3,tmp2,prec);
	  
	}
      
      fclose(infile);
      printf("We processed %lu zeros in %s.\n",zero_count,fname);
      
    }

  fclose(infile1);
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



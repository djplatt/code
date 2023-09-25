/*
Make up a degree 2 L-function by multiplying the quadratic character
mod 3 and something else together.
*/
#include "errno.h"
#include <acb_poly.h>
#include "glfunc.h"
#include "glfunc_internals.h"

#define DIGITS 20
#define MAX_Q 3710369067404ll
#define MAX_FACS 12 // 4*3*5*7...*37 > 3710369067404

long int read_factors(long int *factors, FILE *infile)
{
  long int q;
  if(fscanf(infile,"%ld",&q)!=1)
    return 0;
  int f_ptr=0;
  while(true)
    {
      if(fscanf(infile,"%ld",factors+f_ptr)!=1)
	{
	  printf("Error reading factor file.\n");
	  return 0;
	}
      if(factors[f_ptr]==0)
	break;
      f_ptr++;
    }
  fscanf(infile,"\n");
  return q;
}

// compute the Euler poly for p
// with L the product of primitive char mod 3 and primitive real character mod Q
void lpoly_callback(acb_poly_t poly, uint64_t p, int d __attribute__((unused)), int64_t prec, void *vfactors)
{
  static bool init = false;
  static acb_poly_t poly_1,poly_1_1,poly_1_m1,poly_1_2_1,poly_1_m2_1,poly_1_0_1,poly_1_0_m1;
  static fmpz_t fa,fq;
  if(!init)
    {
      init=true;
      fmpz_init(fa);
      fmpz_init(fq);
      acb_poly_init(poly_1);
      acb_poly_init(poly_1_1);
      acb_poly_init(poly_1_m1);
      acb_poly_init(poly_1_2_1);
      acb_poly_init(poly_1_m2_1);
      acb_poly_init(poly_1_0_m1);
      acb_poly_one(poly_1); // 1
      acb_poly_one(poly_1_1);
      acb_poly_set_coeff_si(poly_1_1,1,1); // 1+T
      acb_poly_one(poly_1_m1);
      acb_poly_set_coeff_si(poly_1_m1,1,-1); // 1-T
      acb_poly_one(poly_1_2_1);
      acb_poly_set_coeff_si(poly_1_2_1,1,2);
      acb_poly_set_coeff_si(poly_1_2_1,2,1); // 1+2T+T^2
      acb_poly_one(poly_1_m2_1);
      acb_poly_set_coeff_si(poly_1_m2_1,1,-2);
      acb_poly_set_coeff_si(poly_1_m2_1,2,1); // 1-2T+T^2
      acb_poly_one(poly_1_0_m1);
      acb_poly_set_coeff_si(poly_1_0_m1,2,-1); // (1-T^2)
    }

  long int *factors=vfactors;
  fmpz_set_si(fa,p);


  int ch=1;
  int fac_ptr=0;
  if(factors[0]==4)
    {
      fac_ptr++;
      if(p==2)
	ch=0;
      else
	{
	  if((p%4)==3)
	    ch=-1;
	}
    }
  if(factors[0]==8)
    {
      fac_ptr++;
      if(p==2)
	ch=0;
      else
	{
	  int res=p%8;
	  if((res==3)||(res==5))
	    ch=-1;
	}
    }
  if(factors[0]==-8)
    {
      fac_ptr++;
      if(p==2)
	ch=0;
      else
	{
	  ch=-1;
	  int res=p%8;
	  if((res==1)||(res==3))
	    ch=1;
	}
    }
  while(factors[fac_ptr])
    {
      if(p==factors[fac_ptr])
	{
	  ch=0;
	  break;
	}
      fmpz_set_si(fq,factors[fac_ptr]);
      ch*=fmpz_jacobi(fa,fq);
      fac_ptr++;
    }
  if(ch==0) // poly = 1
    {
      int p3=p%3;
      if(p3==0)
	acb_poly_set(poly,poly_1);
      else
	{
	  if(p3==1) // 1 mod 3 -> 1-T
	    acb_poly_set(poly,poly_1_m1);
	  else // 2 mod 3 -> 1+T
	    acb_poly_set(poly,poly_1_1);
	}
    }
  else
    {
      if(ch==1) // 1-T
	{
	  int p3=p%3;
	  if(p3==0) // 1
	    acb_poly_set(poly,poly_1_m1); // (1-T) (1) = (1-T)
	  else
	    {
	      if(p3==1) // 1-T
		acb_poly_set(poly,poly_1_m2_1); // (1-T) (1-T) = (1-2T+T^2)
	      else // 1+T
		acb_poly_set(poly,poly_1_0_m1); // (1-T) (1+T) = (1-T^2)
	    }
	}
      else
	// ch==-1 => 1+T
	{
	  int p3=p%3;
	  if(p3==0) // 1
	    acb_poly_set(poly,poly_1_1); // (1+T) (1) = (1+T)
	  else
	    {
	      if(p3==1) // 1-T
		acb_poly_set(poly,poly_1_0_m1); // (1+T) (1-T) = (1-T^2)
	      else // 1+T
		acb_poly_set(poly,poly_1_2_1); // (1+T) (1+T) = (1+2T+T^2)
	    }
	}
    }
  //if(p<20) {printf("ch was %d poly for prime %d is ",ch,p);acb_poly_printd(poly,20);printf("\n");}
  return;
}


int main (int argc, char**argv)
{
  printf("Command Line:- %s",argv[0]);
  for(int i=1;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=2)
    {
      printf("Usage:- %s <infile>\n",argv[0]);
      return 0;
    }

  FILE *infile=fopen(argv[1],"r");
  if(errno)
    {
      char err_str[1024];
      perror("Fatal error in main:-");
      return 0;
    }

  Lfunc_t L;
  double mus[]={1.0,1.0};

  long int factors[MAX_FACS];

  long int Q;
  // we have a degree 2 L-function with alg=anal so normalisation = 0.0
  while(Q=read_factors(factors,infile))
    {
      printf("Doing even character for modulus %ld: ",Q);
      Lerror_t ecode;
      L=Lfunc_init(2,Q*3,0.0,mus,&ecode);
      Lfunc *lf=(Lfunc *) L;
	
      if(fatal_error(ecode))
	{
	  fprintf(stderr,"Fatal error processing modulus %ld\n",Q);
	  fprint_errors(stderr,ecode);
	  printf("\n");
	  Lfunc_clear(L);
	  fflush(stdout);
	  continue;
	}

      lf->self_dual=YES; // product of two self duals is self dual

      long int nmax=Lfunc_nmax(L);
      //printf("nmax=%ld\n",nmax);
      Lfunc_reduce_nmax(L,3*nmax/4);

      ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, factors);

      if(fatal_error(ecode))
	{
	  fprintf(stderr,"Fatal error processing modulus %ld\n",Q);
	  fprint_errors(stderr, ecode);
	  printf("\n");
	  fflush(stdout);
	  Lfunc_clear(L);
	  continue;
	}
      
      // do the computation
      ecode|=Lfunc_compute(L);

      if(fatal_error(ecode))
	{
	  fprintf(stderr,"Fatal error processing modulus %ld\n",Q);
	  fprint_errors(stderr,ecode);
	  printf("\n");
	  fflush(stdout);
	  Lfunc_clear(L);
	  continue;
	}

      arb_t *zeros=(arb_t *)Lfunc_zeros(L,0);
      printf("First zero at:- ");
      arb_printd(zeros[0],35);
      printf("\n");
      fflush(stdout);

      Lfunc_clear(L);

      // print any warnings collected along the way
      if(ecode!=ERR_SUCCESS)
	{
	  fprintf(stderr,"Non-fatal error processing modulus %ld\n",Q);
	  fprint_errors(stderr,ecode);
	}
    }
  fclose(infile);

  return 0;
}



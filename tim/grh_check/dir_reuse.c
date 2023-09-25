/*
Make up a degree 2 L-function by multiplying the quadratic character
mod 3 and something else together.
*/
#include "errno.h"
#include <acb_poly.h>
#include "glfunc.h"
#include "glfunc_internals.h"

#define DIGITS 20
#define MAX_Q 300000000
#define MAX_FACS 8 // 4*3*5*7...*23 > 300000000

int read_factors(int *factors, FILE *infile)
{
  int q;
  if(fscanf(infile,"%d",&q)!=1)
    return 0;
  int f_ptr=0;
  while(true)
    {
      if(fscanf(infile,"%d",factors+f_ptr)!=1)
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
// with L the product of non-principal characters mod 5 and 7
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

  int *factors=vfactors;
  fmpz_set_si(fa,p);


  acb_poly_t pq;
  int ch=1;
  acb_poly_init(pq);
  acb_poly_one(pq);
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
  double mus[]={0.0,1.0};

  int factors[MAX_FACS];

  int Q,first_Q;
  bool first=true;

  // we have a degree 2 L-function with alg=anal so normalisation = 0.0
  while(Q=read_factors(factors,infile))
    {
      Lerror_t ecode;
      Lfunc *lf;
      
      if(first)
	{
	  first=false;
	  first_Q=Q;
	  L=Lfunc_init(2,Q*3,0.0,mus,&ecode);
	  if(fatal_error(ecode))
	    {
	      fprint_errors(stderr,ecode);
	      return 0;
	    }

	  lf=(Lfunc *) L;
	  lf->self_dual=YES;
	}
      else
	{
	  if(Q>first_Q)
	    // we need Q to be decreasing because the error bounds
	    // improve with lower Q (so no need to re-compute them)
	    {
	      printf("Fatal error: We expect the first Q seen to be the largest. Exiting.\n");
	      return 0;
	    }
	  lf->conductor=3*Q; // conductor changes with Q
	  lf->nmax_called=false; // sorts out re-init
	  //arb_zero(lf->buthe_Wf); // this is computed additively per a(n)
	  //for(int n=0;n<lf->M;n++)
	  //acb_set_ui(lf->ans[n],1); // the a(n) are computed mult.
	}
      
      printf("Doing even character for modulus %d\n",Q);
      /*
      printf("Normalisation:- %f\n",lf->normalisation);
	      
      for(int n=0;n<10;n++)
	{
	  printf("a(%d)=",n+1);acb_printd(lf->ans[n],20);printf("\n");
	}
      */

      ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, factors);
      if(fatal_error(ecode))
	{
	  fprint_errors(stderr, ecode);
	  return 0;
	}
      /*
      printf("Normalisation:- %f\n",lf->normalisation);

      for(int n=0;n<10;n++)
	{
	  printf("a(%d)=",n+1);acb_printd(lf->ans[n],20);printf("\n");
	}
      */
      
      // do the computation
      ecode|=Lfunc_compute(L);
      if(fatal_error(ecode))
	{
	  fprint_errors(stderr,ecode);
	  return 0;
	}

      /*      
      printf("First 10 zeros\n");
      // we could use Lfunc_zeros(L, 1) for the dual L-function
      arb_srcptr zeros=Lfunc_zeros(L, 0);
      for(int i  = 0; i < 10; ++i) {
	printf("Zero %d = ", i);
	arb_printd(zeros+i, DIGITS);
	printf("\n");
      }

      Lplot_t *Lp=Lfunc_plot_data(L,0,10.0,21);
      double t=0.0;
      for(int n=0;n<Lp->n_points;n++,t+=Lp->spacing)
	printf("%f %f\n",t,Lp->points[n]);
      */
 
      // print any warnings collected along the way
      fprint_errors(stderr,ecode);
    }

  Lfunc_clear(L);
  fclose(infile);
  return 0;
}



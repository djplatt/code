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
      int p7=p%7;
      if(p7==0)
	acb_poly_set(poly,poly_1);
      else
	{
	  if((p7==1)||(p7==2)||(p7==4)) // 1 mod 3 -> 1-T
	    acb_poly_set(poly,poly_1_m1);
	  else // 2 mod 3 -> 1+T
	    acb_poly_set(poly,poly_1_1);
	}
    }
  else
    {
      if(ch==1) // 1-T
	{
	  int p7=p%7;
	  if(p7==0) // 1
	    acb_poly_set(poly,poly_1_m1); // (1-T) (1) = (1-T)
	  else
	    {
	      if((p7==1)||(p7==2)||(p7==4)) // 1-T
		acb_poly_set(poly,poly_1_m2_1); // (1-T) (1-T) = (1-2T+T^2)
	      else // 1+T
		acb_poly_set(poly,poly_1_0_m1); // (1-T) (1+T) = (1-T^2)
	    }
	}
      else
	// ch==-1 => 1+T
	{
	  int p7=p%7;
	  if(p7==0) // 1
	    acb_poly_set(poly,poly_1_1); // (1+T) (1) = (1+T)
	  else
	    {
	      if((p7==1)||(p7==2)||(p7==4)) // 1-T
		acb_poly_set(poly,poly_1_0_m1); // (1+T) (1-T) = (1-T^2)
	      else // 1+T
		acb_poly_set(poly,poly_1_2_1); // (1+T) (1+T) = (1+2T+T^2)
	    }
	}
    }
  //if(p<20) {printf("ch was %d poly for prime %d is ",ch,p);acb_poly_printd(poly,20);printf("\n");}
  return;
}

double normalised(Lfunc *,uint64_t, uint64_t,double);

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

  int Q;
  // we have a degree 2 L-function with alg=anal so normalisation = 0.0
  while(Q=read_factors(factors,infile))
    {
      printf("Doing even character for modulus %d: ",Q);
      Lerror_t ecode;
      L=Lfunc_init(2,Q*7,0.0,mus,&ecode);
      Lfunc *lf=(Lfunc *) L;
	
      if(fatal_error(ecode))
	{
	  fprintf(stderr,"Fatal error processing modulus %d\n",Q);
	  fprint_errors(stderr,ecode);
	  printf("\n");
	  Lfunc_clear(L);
	  continue;
	}

      lf->self_dual=YES; // product of two self duals is self dual

      ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, factors);

      if(fatal_error(ecode))
	{
	  fprintf(stderr,"Fatal error processing modulus %d\n",Q);
	  fprint_errors(stderr, ecode);
	  printf("\n");
	  Lfunc_clear(L);
	  continue;
	}
      
      // do the computation
      ecode|=Lfunc_compute(L);

      if(fatal_error(ecode))
	{
	  fprintf(stderr,"Fatal error processing modulus %d\n",Q);
	  fprint_errors(stderr,ecode);
	  printf("\n");
	  Lfunc_clear(L);
	  continue;
	}

      arb_t *zeros=(arb_t *)Lfunc_zeros(L,0);
      printf("First zero at:- ");
      arb_printd(zeros[0],35);
      printf("\n");

      Lfunc_clear(L);

      // print any warnings collected along the way
      if(ecode!=ERR_SUCCESS)
	{
	  fprintf(stderr,"Non-fatal error processing modulus %d\n",Q);
	  fprint_errors(stderr,ecode);
	}
    }
  fclose(infile);

  return 0;
}


/*

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
  //double mus[]={0.0,1.0};

  int factors[MAX_FACS];

  int Q;


  // we have a degree 2 L-function with alg=anal so normalisation = 0.0
  while(Q=read_factors(factors,infile))
    {
      printf("Doing even character for modulus %d: ",Q);
      Lerror_t ecode;
      Lparams_t Lp;
      Lp.degree=2;
      Lp.conductor=7*Q;
      Lp.normalisation=0.0;
      Lp.mus=(double *)malloc(sizeof(double)*2);
      Lp.mus[0]=0.0;
      Lp.mus[1]=1.0;
      Lp.target_prec = DEFAULT_TARGET_PREC;
      Lp.rank = DK;
      Lp.self_dual = YES;
      Lp.cache_dir = ".";
      Lp.gprec = 0; // We will try to do something sensible
      Lp.wprec = 0; // ditto
      L=Lfunc_init_advanced(&Lp,&ecode);
      //L=Lfunc_init(2,Q*3,0.0,mus,&ecode);
      Lfunc *lf=(Lfunc *) L;
	
      if(fatal_error(ecode))
	{
	  fprintf(stderr,"Fatal error processing modulus %d\n",Q);
	  fprint_errors(stderr,ecode);
	  return 0;
	}

      ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, factors);

      if(fatal_error(ecode))
	{
	  fprintf(stderr,"Fatal error processing modulus %d\n",Q);
	  fprint_errors(stderr, ecode);
	  return 0;
	}
      
      // do the computation
      ecode|=Lfunc_compute(L);

      if(fatal_error(ecode))
	{
	  fprintf(stderr,"Fatal error processing modulus %d\n",Q);
	  fprint_errors(stderr,ecode);
	  return 0;
	}
      
      arb_t *zeros=(arb_t *)Lfunc_zeros(L,0);
      printf("First zero at:- ");
      arb_printd(zeros[0],35);
      printf("\n");

      
      double zeros_7[]={4.47573828372868313197462848719, 6.84549171249137726783979779478, 11.16018454311952965510181826588, 12.48960334303313423775824003050, 15.11288225874376865962836973665, 16.802876475728853982349899698677, 19.61187805669003102874335093663, 21.89991370331832248627181534368, 23.16297179973657863840493856183, 24.49884755535632897308961801424, 27.36147751900128252344458972943, 29.17887966636347069628776553811, 30.72189964454757558763990410272, 32.4758875915420758078446454007, 34.11124078303733952175233551521, 35.311847458977969626587230969278, 38.36463733603136596844459180378, 39.48257874666192172268872121810, 40.7103286516539188858032460596, 42.30485758764265058072877014058, 43.93185428742714565723041221515, 46.19432742740581606201873969732, 47.29858534027669804962753511398, 0.0}; // zeros of L_7
      arb_t tmp,tmp1,small;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(small);
      arb_set_d(small,1e-10);
      arb_set_d(tmp,zeros_7[0]); // first zero of L_3
      for(int three_ptr=1,i=0,j=0;;i++)
	{
	  if(arb_is_zero(zeros[i])) // end of zeros
	    break;
	  arb_sub(tmp1,zeros[i],tmp,100);
	  arb_abs(tmp1,tmp1);
	  arb_sub(tmp1,tmp1,small,100);
	  if(arb_is_negative(tmp1)) // its (very close to) a zero of L_3
	    arb_set_d(tmp,zeros_7[three_ptr++]);
	  else
	    {
	      printf("zero %d at:- ",j++);
	      arb_printd(zeros[i],35);
	      printf("\n");
	    }
	}
      arb_clear(tmp);
      arb_clear(tmp1);
      arb_clear(small);
      
      
      for(int n=0;n<=lf->fft_NN/OUTPUT_RATIO+lf->fft_NN/TURING_RATIO;n++)
	printf("%f %g\n",n/256.0,normalised(lf,0,n,n/256.0));
      
      arb_t sum;
      arb_init(sum);
      arb_add(sum,lf->buthe_Wf,lf->buthe_Winf,25);
      arb_sub(sum,sum,lf->buthe_Ws,25);
      printf("Ws = ");arb_printd(lf->buthe_Ws,20);
      printf("\nWinf = ");arb_printd(lf->buthe_Winf,20);
      printf("\nWf = ");arb_printd(lf->buthe_Wf,20);
      printf("\nS = ");arb_printd(sum,20);printf("\n");
      arb_clear(sum);

      

      Lfunc_clear(L);

      // print any warnings collected along the way
      if(ecode!=ERR_SUCCESS)
	{
	  fprintf(stderr,"Non-fatal error processing modulus %d\n",Q);
	  fprint_errors(stderr,ecode);
	}
    }
  fclose(infile);

  return 0;
}

*/

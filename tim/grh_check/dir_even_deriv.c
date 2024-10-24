/*
Make up a degree 2 L-function by multiplying the quadratic character
mod 3 and an even character mod q together.
Compute Lambda'/Lambda(1) and compare with R log q
*/
#include "errno.h"
#include <flint/acb_poly.h>
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
  //fscanf(infile,"\n");
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

  long int *factors=(long int *)vfactors;
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

  int64_t prec=200;
  printf("Command Line:- %s",argv[0]);
  for(int i=1;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=2)
    {
      printf("Usage:- %s <q file>\n",argv[0]);
      return 0;
    }

  FILE *infile=fopen(argv[1],"r");
  if(errno)
    {
      char err_str[1024];
      perror("Fatal error in main:-");
      return 0;
    }
  /*
  FILE *threefile=fopen(argv[1],"r");
  if(errno)
    {
      char err_str[1024];
      perror("Fatal error in main:-");
      return 0;
    }
  */
  
  arb_t R;
  arb_init(R);
  arb_set_ui(R,9645908801);
  arb_div_ui(R,R,1000000000,prec);
  arb_t RlogQ;
  arb_init(RlogQ);
  
  arb_t Lam3,Lam3_dash;
  arb_init(Lam3);arb_init(Lam3_dash);
  // these were computed off line by glfunc
  arb_load_str(Lam3,"40d6dd2411b11d7eb6b55f92ddaae4846ac9b83 -9c a49a929 -9b");
  arb_load_str(Lam3_dash,"eaefa4ab8776b96a29456fec835f74ca640cf -9a f27ed1b -9a");

  printf("Lam_3(1) = ");arb_printd(Lam3,20);
  printf("\nLam'_3(1) = ");arb_printd(Lam3_dash,20);printf("\n");
  
  double mus[]={0.0,1.0}; // L3 is odd, other is even

  long int factors[MAX_FACS];

  long int Q;

  acb_t Lam,Lam_dash,Lamq,Lamq_dash;
  acb_init(Lam);
  acb_init(Lam_dash);
  acb_init(Lamq);
  acb_init(Lamq_dash);

  acb_t ctmp;
  acb_init(ctmp);
  arb_t rtmp;
  arb_init(rtmp);
  
  // we have a degree 2 L-function with alg=anal so normalisation = 0.0
  while(Q=read_factors(factors,infile))
    {
      //printf("Doing even character for modulus %ld: ",Q);
      Lerror_t ecode;
      Lfunc_t L;
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

      ecode|=Lfunc_special_value_choice(Lam, Lam_dash, L, 1, 0.0,true,true);
      if(fatal_error(ecode))
	{
	  fprint_errors(stderr,ecode);
	  return 0;
	}

      acb_div_arb(Lamq,Lam,Lam3,prec);
      acb_mul_arb(ctmp,Lamq,Lam3_dash,prec);
      acb_sub(ctmp,Lam_dash,ctmp,prec);
      acb_div_arb(Lamq_dash,ctmp,Lam3,prec);

      
      //printf("Lam(1) = ");acb_printd(Lamq, DIGITS);printf("\n");
      //printf("Lam'(1) = ");acb_printd(Lamq_dash, DIGITS);printf("\n");
      acb_div(ctmp,Lamq_dash,Lamq,prec);
      printf("Q : %lu : Lam'/Lam(1) = ",Q);acb_printd(ctmp,DIGITS);
      arb_log_ui(RlogQ,Q,prec);
      arb_mul(RlogQ,RlogQ,R,prec);
      arb_sub(rtmp,acb_realref(ctmp),RlogQ,prec);
      if(arb_is_negative(rtmp))
	printf(" passed\n");
      else
	printf(" failed\n");

      Lfunc_clear(L);
  
      // print any warnings collected along the way
      if(ecode!=ERR_SUCCESS)
	{
	  fprintf(stderr,"Non-fatal error processing modulus %ld\n",Q);
	  fprint_errors(stderr,ecode);
	}
      fflush(stdout);
      
    }     

  acb_clear(Lam);acb_clear(Lam_dash);
  acb_clear(Lamq);acb_clear(Lamq_dash);
  arb_clear(Lam3);arb_clear(Lam3_dash);
  acb_clear(ctmp);arb_clear(R);arb_clear(RlogQ);
  arb_clear(rtmp);
  
  fclose(infile);
  
  return 0;
}



// Compute factoriastion data for G3 curves
// we want 1/pi(x) sum_{p\leq x} a_p^2/p
#include "stdio.h"
#include "inttypes.h"
#include "string.h"
#include "flint/acb_poly.h"
#include "smalljac.h"
#include "../andy_elliptic/pari.c"
//#include "../generic/defines.h"

#define g2_ell_arb_prec (100)
#define MAX_EULER_FACTOR (100)

//**********************************************************************
// utility routines by Andy Booker
//**********************************************************************


// structure to hold the bad primes data
#define MAX_DEGREE (6) // Genus 3 curves
typedef struct {
	long int p,f[MAX_DEGREE+1];
} bad_lfactors_t;

// there can't be more bad primes than this
bad_lfactors_t bad_lfactors[PRIMORIAL_MAX_COND+1];

extern "C" int i_poly_parse(long *,int,char *);

//double sum;
//uint64_t pi_x;
// handle an Euler Polynomial of degree DEGREE
// polynomial coefficients are long ints
int do_lpoly_smalljac(smalljac_Qcurve_t C,unsigned long p,
	     int good,long aa[],int n,void *arg,int DEGREE) 
{
  if(DEGREE>MAX_DEGREE)
    {
      fprintf(stderr,"DEGREE exceeded MAX_DEGREE in do_lpoly. Exiting.\n");
      exit(0);
    }

  if (good) 
    {
      //sum+=aa[0]*aa[0]/(double) p;
      //pi_x++;
      //printf("prime: %lu sum now %10.8e\n",p,sum);
      return(0);
    }
  else 
    {
      uint64_t i;
      for (i=0;bad_lfactors[i].p;i++)
	if (bad_lfactors[i].p == p)
	  {
	    //sum+=bad_lfactors[i].f[1]*bad_lfactors[i].f[1]/(double) p;
	    //pi_x++;
	    //printf("bad prime: %lu sum now %10.8e\n",p,sum);
	    break;
	  }
      if (!bad_lfactors[i].p) {
	fprintf(stderr,"unexpected bad reduction at p=%d. Exiting.\n",p);
	exit(1);
      }
    }
  return 1;
}


int do_lpoly_g3(smalljac_Qcurve_t C,unsigned long p,
		int good,long aa[],int n,void *arg) 
{
  printf("In do_lpoly_g3 with p=%lu\n",p);
  return do_lpoly_smalljac(C,p,good,aa,n,arg,6);
}

bool g3_parse_line(char *str)
{
  smalljac_Qcurve_t curve;
  char *hash,*s,*pq;
  int64_t disc,cond;
  int k,p,d,i;
  if (!(s=strtok(str,":")) || sscanf(s,"%ld",&disc) != 1) 
    {
      fprintf(stderr,"Bad discriminant.\n");
      return false;
    }
  if (!(s=strtok(NULL,":")) || sscanf(s,"%ld",&cond)!=1)
    {
      fprintf(stderr,"Bad conductor.\n");
      return false;
    }

  if(!(hash=strtok(NULL,":"))) 
    {
      fprintf(stderr,"Bad hash.\n");
      return false;
    }
  if (!(pq=strtok(NULL,":"))) return false; // string describing curve to smalljac
  fprintf(stderr,"Calling smalljac_init with %s\n",pq);fflush(stderr);
  if ( !(curve=smalljac_curve_init(pq,&k)) )
    {
      fprintf(stderr,"smalljac_curve_init failed.\n");
      return false;
    }
  if(!strtok(NULL,":")) return false; // not sure what this field is
  if (!(pq=strtok(NULL," "))) return false;
  fprintf(stderr,"Bad prime data = %s\n",pq);fflush(stderr);
  if(pq[1]=='0')
    {
      fprintf(stderr,"Don't know what to do with prime 2 for this curve.\n");
      return false;
    }
  p = 2, k = 0;
  for (s=strtok(pq+1,",");s;s=strtok(NULL,","))
    {
      for (;p<=disc;p++)
	if (disc % p == 0) {
	  do disc /= p; while (disc % p == 0);
	  bad_lfactors[k].p = p;
	  d = i_poly_parse(bad_lfactors[k].f,6,s);
	  for (i=d+1;i<=6;i++)
	    bad_lfactors[k].f[i] = 0;
	  k++;
	  break;
	}
    }
  bad_lfactors[k].p = 0;
  if (disc != 1) return false;

  //sum=0.0;
  //pi_x=0;
  smalljac_Lpolys(curve,1,MAX_EULER_FACTOR,0,do_lpoly_g3,NULL);
  smalljac_Qcurve_clear(curve);
  //printf("Normalised sum =%10.8e\n",sum/(double) pi_x);
  return true;
}


int main(int argc, char **argv)
{
  fprintf(stderr,"Command line:- %s",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    fprintf(stderr," %s",argv[i]);
  fprintf(stderr,"\n");
  if(argc!=2)
    {
      fprintf(stderr,"Usage:- %s <g3 string>.\n",argv[0]);
      return 0;
    }

  g3_parse_line(argv[1]);

  return 0;
}

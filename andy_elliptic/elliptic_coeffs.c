//
// output (analytic, i.e. *n^{-1/2}) Diriclet coeefficients for elliptic curve
// described by argv[1]
//
// feed output into code/glfunc/artin and data/modforms/gdata/mf/mf.2
//
// DJ Platt Jan 2017
//
// cribbed from Andy Booker's plot2.c code
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/file.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <math.h>
#include <acb.h>
#include "smalljac.h"
#include "pari.c"
#include "hashcoeffs.c"

extern int i_poly_parse(long *,int,char *);
//extern acb_t *fft_init(long,long);
//extern void fft(acb_t *);
//extern void ifft(acb_t *);
//extern void K0(arb_t,arb_srcptr);

#define prec 300
#define degree 2
#if degree == 2
#define A 4
#define stride 16
#else // degree = 4
#define A 5
#define stride 15
#endif
#define sample_points 1024
#define upsample_factor 16
#define big_sample_points (upsample_factor*sample_points)
#define bigA (A*upsample_factor)
#define fpoints (22*bigA)
#define plot_delta (double)stride/bigA


static int *a,*ftable;
static smalljac_curve_t curve;
static long maxn;
static acb_t *data,*datacopy,*bigdata;
static int epsilon,Rank;
static long hash;
static arb_t scale;
static arb_t *zero;
static int nzeros,mzeros;
static arb_t f[fpoints];

/* invert the polynomial p to precision k */
static void inverse(int *ans,int *p,int k) {
	int i,j,c[32];

	c[0] = 1;
	for (j=1;j<=k;j++) c[j] = 0;
	for (i=0;i<=k;i++) {
		ans[i] = c[0];
		for (j=1;j<=k-i;j++)
			c[j-1] = c[j]-p[j]*ans[i];
	}
}

// floor(log(x)/log(p))
static inline int logp(unsigned long x,unsigned long p) {
	int k;
	for (k=0;x>=p;k++)
		x /= p;
	return k;
}

static struct {
	long p,f[degree+1];
} bad_lfactors[8];

static int recordlpoly(smalljac_Qcurve_t C,unsigned long p,
int good,long aa[],int n,void *arg) {
	int f[64],c[64],i,j,k;
	unsigned long q;

	f[0] = 1; for (j=1;j<64;j++) f[j] = 0;
	if (good) {
#if degree == 2
		f[1] = aa[0];
		f[2] = p;
#else // degree 4
		f[1] = aa[0];
		f[2] = aa[1];
		f[3] = f[1]*p;
		f[4] = p*p;
#endif
	} else {
		for (i=0;bad_lfactors[i].p;i++)
			if (bad_lfactors[i].p == p) {
				for (j=0;j<=degree;j++)
					f[j] = bad_lfactors[i].f[j];
				break;
			}
		if (!bad_lfactors[i].p) {
			fprintf(stderr,"unexpected bad reduction\n");
			exit(1);
		}
	}

	k = logp(maxn,p);
	inverse(c,f,k);
	for (j=1,q=1;j<=k;j++)
		q *= p, a[q] = c[j];
	return 1;
}

static void allocate(long n) {
	int k;
	static long nallocated;

	if (n <= nallocated) return;

	// allocate space for coefficients and factor table
	if (nallocated) a++, ftable++;
	a = (int *)realloc(a,n*sizeof(*a))-1;
	ftable = (int *)realloc(ftable,n*sizeof(*ftable))-1;
	nallocated = n;

  // sieve to compute first factor table */
	for (n=1;n<=nallocated;n++)
		ftable[n] = n;
	for (k=(int)sqrtl((long double)nallocated);k>=2;k--)
		for (n=k+k;n<=nallocated;n+=k)
			ftable[n] = k;
}

static long computehash(void) {
	int i;
	long hash=0;

	for (i=0;i<hashlength;i++)
		hash = (long)((hash+(__int128_t)hashcoeff[i].c
		* a[hashcoeff[i].p]) % hashmodulus);
	if (hash < 0) hash += hashmodulus;
	return hash;
}

void output_arf(arf_t x)
{
  static int init=0;
  static fmpz_t m,e;
  if(!init)
    {
      init=1;
      fmpz_init(m);
      fmpz_init(e);
    }
  arf_get_fmpz_2exp(m,e,x);
  fmpz_print(m);
  printf(" ");
  fmpz_print(e);
}
  

void output_arb(arb_t x)
{
  static int init=0;
  static arb_t tmp;
  static arf_t tmp1;
  if(!init)
    {
      init=1;
      arb_init(tmp);
      arf_init(tmp1);
    }
  arb_get_mid_arb(tmp,x);
  arb_get_ubound_arf(tmp1,tmp,prec);
  //printf("outputting centre = ");arf_printd(tmp1,10);printf("\n");
  output_arf(tmp1);
  printf(" ");
  arb_get_rad_arb(tmp,x);
  arb_get_ubound_arf(tmp1,tmp,prec);
  output_arf(tmp1);
}

void output_acb(acb_t z)
{
  output_arb(acb_realref(z));
  printf(" ");
  output_arb(acb_imagref(z));
  printf("\n");
}

void output_dirichlet_coefficient(int an)
{
  static int init=0;
  static acb_t z;
  if(!init)
    {
      init=1;
      acb_init(z);
      arb_zero(acb_imagref(z));
    }
  arb_set_si(acb_realref(z),an);
  output_acb(z);
}

void output_dirichlet_coefficient_arb(arb_t an)
{
  static int init=0;
  static acb_t z;
  if(!init)
    {
      init=1;
      acb_init(z);
      arb_zero(acb_imagref(z));
    }
  arb_set(acb_realref(z),an);
  output_acb(z);
}



// truncate S at cutoff*x
#define cutoff 33

void compute_coeffs(double xmax) {
  static int init=0;
  static arb_t tmp,tmp1,tmp2;
  if(!init)
    {
      init=1;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
    }
	long p,n,q;

	maxn = (long)(cutoff*xmax);
	// compute at least up to 8192 so we can work out the hash
	if (maxn < 8192) maxn = 8192;
	allocate(maxn);
	a[1] = 1;
	output_dirichlet_coefficient(1);
	smalljac_Lpolys(curve,1,maxn,0,recordlpoly,(void *)0);
	for (n=2;n<=maxn;n++) {
		p = ftable[n];
		for (q=p;n/q%p==0;q*=p);
		a[n] = a[q] * a[n/q];
		arb_sqrt_ui(tmp,n,prec);
		arb_set_si(tmp1,a[n]);
		arb_div(tmp2,tmp1,tmp,prec);
		output_dirichlet_coefficient_arb(tmp2);
	}
	hash = computehash();
}

int main(int argc,char *argv[]) 
{
  int i,k,d,line,kprocs=0;
  char buf[256],lpolys[256],*PQ,*s,last_isog[8];
  long disc,cond,p,last_cond=0;
  FILE *infile;
  arb_t c,t;
  acb_t z;  

  if (argc!=2)
    {
      printf("Usage:- %s <curve definition>.\n",argv[0]);
      return(0);
    }
  
  allocate(cutoff*600); // conductors up to 360000
  if (!(s=strtok(argv[1]," ")) || sscanf(s,"%ld",&cond) != 1) return(0);
  if (!(s=strtok(NULL," ")) ) return(0);
  if (cond == last_cond && !strcmp(s,last_isog)) return(0);
  last_cond = cond; strcpy(last_isog,s);
  if (!(s=strtok(NULL," ")) ) return(0);
  if (!(s=strtok(NULL," ")) ) return(0);
  PQ = s;
  if (!(s=strtok(NULL," ")) || sscanf(s,"%d",&Rank) != 1) return(0);
  epsilon = (Rank & 1) ? -1 : 1;
  disc = cond;

  //printf("Rank = %d\n",Rank);
  //printf("calling smalljac with %s\n",PQ);exit(0);

  if ( !(curve=smalljac_curve_init(PQ,&k)) )
    return(0);
  p = 2, k = 0;
  for (;p<=disc;p++)
    if (disc % p == 0) {
      do disc /= p; while (disc % p == 0);
      bad_lfactors[k].p = p;
      bad_lfactors[k].f[0] = 1;
      sprintf(lpolys,"ellap(ellinit(%s),%ld)",PQ,p);
      bad_lfactors[k].f[1] = -parilong(lpolys);
      bad_lfactors[k].f[2] = 0;
      k++;
    }
  bad_lfactors[k].p = 0;
  if (disc != 1) {
    fprintf(stderr,"something has gone wrong with bad L-factors\n");
    return 1;
  }
  compute_coeffs(sqrt((double)cond));
  smalljac_Qcurve_clear(curve);
  return 0;
}
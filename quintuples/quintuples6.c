/*
 *
 * Count how many triples of the form {a,b,d} with a<b<d, b^6<=d<=b^7.7
 * can be extended to {a,b,c,d} with b<c<d
 * We expect the answer to be zero
 */
#include "inttypes.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"

inline uint64_t gcd (uint64_t a, uint64_t b)
/* Euclid algorithm gcd */
{
	uint64_t c;
	while(a!=0)
	{
		c=a;
		a=b%a;
		b=c;
	};
	return(b);
};


#define MAX_B (1300000000L)
#define MAX_FACS (17) // 2*3*..59>49e18

// structure to hold factor information
typedef struct
{
  //uint64_t pr;   // = 0 if no primitive root
  //uint64_t phi;  // Euler phi(q)
	uint64_t num_facs;  // number of prime power factors
	uint64_t primes[MAX_FACS]; // the prime factors p of n
	uint64_t powers[MAX_FACS];  // the powers of p
	uint64_t facs[MAX_FACS]; // the p^power
} factor_t;

// read factorisations from factors.gp
int old_read_factorisation(FILE *infile, uint64_t *r, factor_t *f)
{
  if(fscanf(infile,"%lu %lu",r,&(f->num_facs))!=2)
    return(0);
  uint64_t ptr;
  for(ptr=0;ptr<f->num_facs;ptr++)
    {
      if(fscanf(infile,"%lu %lu",f->primes+ptr,f->powers+ptr)!=2)
	{
	  printf("Failed to parse factorisation of %lu. Exiting.\n",r[0]);
	  exit(0);
	}
      uint64_t i;
      f->facs[ptr]=f->primes[ptr];
      for(i=1;i<f->powers[ptr];i++)
	f->facs[ptr]*=f->primes[ptr];
    }
  fscanf(infile,"\n");
  return(1);
}

// read factorisations from factors.gp

int read_factorisation(FILE *infile, uint64_t *r, factor_t *f)
{
  if(fscanf(infile,"%lu %lu",r,&(f->num_facs))!=2)
    return(0);
  uint64_t fptr,pptr,last_p=0,p,pows,ffacs=f->num_facs;
  for(fptr=0,pptr=0;fptr<ffacs;fptr++)
    {
      if(fscanf(infile,"%lu %lu",&p,&pows)!=2)
	{
	  printf("Failed to parse factorisation of %lu. Exiting.\n",r[0]);
	  exit(0);
	}

      if(p==last_p)
	{
	  f->powers[pptr-1]+=pows;
	  uint64_t i;
	  for(i=0;i<pows;i++)
	    f->facs[pptr-1]*=p;
	  f->num_facs--;
	}
      else
	{
	  f->primes[pptr]=p;
	  f->powers[pptr]=pows;
	  f->facs[pptr]=p;
	  uint64_t i;
	  for(i=1;i<pows;i++)
	    f->facs[pptr]*=p;
	  last_p=p;pptr++;
	}
    }
  fscanf(infile,"\n");
  return(1);
}
      

void print_factors(factor_t *f, uint64_t n)
{
  printf("%lu factors as\n",n);
  uint64_t i=0;
  for(i=0;i<f->num_facs;i++)
    printf("  %lu^%lu\n",f->primes[i],f->powers[i]);
}

mpz_t tmp1,tmp2,tmp3,ab1,r,t,s,d,upb,lwb;

mpz_t lwbc,upbc,tc,sc,rc,c;


uint64_t abcalls=0;

// c <-- (x^2-1)/a
// return 1 iff a divides x^2-1
int set_c(mpz_t x, uint64_t a, mpz_t c, mpz_t lim, uint64_t b1)
{
  mpz_div_ui(tmp1,x,b1);
  mpz_mul(tmp1,tmp1,tmp1);
  mpz_sub_ui(tmp1,tmp1,1);
  if(mpz_fdiv_q_ui(c,tmp1,a)==0)
    return(1);
  if(mpz_cmp(c,lim)>=0)
    {
      mpz_set_ui(c,0);
      return(-1);
    }
  mpz_set_ui(c,0);
  return(0);
}

void mpfi_floor(mpz_t res, mpfi_t x, mpfr_ptr rtmp, mpz_t ztmp)
{
  mpfi_get_left(rtmp,x);
  mpfr_get_z(res,rtmp,GMP_RNDD);
  mpfi_get_left(rtmp,x);
  mpfr_get_z(ztmp,rtmp,GMP_RNDD);
  if(mpz_cmp(res,ztmp)!=0)
    {
      printf("Error flooring. Use more precision. Exiting.\n");
      exit(0);
    }
}

uint64_t my_mpz_get_ui(mpz_t x)
{
  if(mpz_cmp_ui(x,0xFFFFFFFFFFFFFFFFLL)>0)
    {
      printf("mpz was too large to fit into ui. Exiting.\n");
      exit(0);
    }
  return(mpz_get_ui(x));
}


mpz_t _P,_Q,_a,_G,_G1,_G2,_B,_B1,_B2,pm1res[2];
int PQm1(mpz_t *res, uint64_t D, mpfi_t sqrt_D, mpfi_t *itmps, mpfr_ptr rtmp, mpz_t ztmp)
{
  mpz_set_ui(_P,0);mpz_set_ui(_Q,1);
  mpfi_floor(_a,sqrt_D,rtmp,ztmp);
  uint64_t i=0;
  mpz_set_ui(_G2,1);
  mpz_set(_G1,_a);
  mpz_set_ui(_B2,0);
  mpz_set_ui(_B1,1);
  while((i==0)||(mpz_cmp_ui(_Q,1)!=0))
    {
      mpz_mul(ztmp,_a,_Q);
      mpz_sub(_P,ztmp,_P);
      mpz_mul(ztmp,_P,_P);
      mpz_neg(ztmp,ztmp);
      mpz_add_ui(ztmp,ztmp,D);
      mpz_div(_Q,ztmp,_Q);
      mpfi_set_z(itmps[0],_P);
      mpfi_add(itmps[1],itmps[0],sqrt_D);
      mpfi_div_z(itmps[0],itmps[1],_Q);
      mpfi_floor(_a,itmps[0],rtmp,ztmp);
      mpz_mul(_G,_a,_G1);
      mpz_add(_G,_G,_G2);
      mpz_set(_G2,_G1);mpz_set(_G1,_G);
      mpz_mul(_B,_a,_B1);
      mpz_add(_B,_B,_B2);
      mpz_set(_B2,_B1);mpz_set(_B1,_B);
      i++;
      /*
      printf("%lu ",i);
      mpz_out_str(NULL,10,_P);printf(" ");
      mpz_out_str(NULL,10,_Q);printf(" ");
      mpz_out_str(NULL,10,_a);printf(" ");
      mpz_out_str(NULL,10,_B);printf(" ");
      mpz_out_str(NULL,10,_G);printf(" ");
      printf("\n");
      */
    }
  mpz_set(res[0],_G2);
  mpz_set(res[1],_B2);
  if(i&1)
    return(-1);
  else      
    return(1);
}
 

#define C_LIM (100)
mpz_t cs[C_LIM];
#define D_LIM (100)
mpz_t ds[D_LIM];
uint64_t cptr,dptr;
#define ZTMPS (4)

/*
void next_xy(mpz_t s, mpz_t t, uint64_t a, uint64_t b, uint64_t ur)
{
  mpz_mul_ui(tmp1,t,ur);
  mpz_mul_ui(tmp2,s,b);
  mpz_add(tmp3,tmp1,tmp2); // t*r+s*b
  mpz_mul_ui(tmp1,t,a);
  mpz_mul_ui(tmp2,s,ur);
  mpz_add(s,tmp1,tmp2); // t*a+r*s
  mpz_set(t,tmp3);
  //printf("x,y set to ");mpz_out_str(NULL,10,s);printf(",");
  //mpz_out_str(NULL,10,t);printf("\n");
}



uint64_t try_it(uint64_t a, uint64_t b, uint64_t ur, int64_t x, int64_t y)
{
  printf("In try_it with %lu %lu %lu %ld %ld\n",a,b,ur,x,y);
  uint64_t count=0;
  mpz_ui_pow_ui(lwb,b,5);
  mpz_ui_pow_ui(upb,b,8); // quicker than b^7.7 ?
  mpz_set_si(s,x);
  mpz_set_si(t,y);
  uint64_t cptr=0;
  set_c(s,a,cs[0]);
  while(mpz_cmp(cs[cptr],lwb)<0)
    {
      printf("found candidate c=");
      mpz_out_str(NULL,10,cs[cptr]);
      printf("\n");
      next_xy(s,t,a,b,ur);
      cptr++;
      if(cptr==C_LIM)
	{
	  printf("Fatal error. Must increase C_LIM. Exiting.\n");
	  exit(0);
	}
      set_c(s,a,cs[cptr]);
    }
  mpz_set(d,cs[cptr]);
  while(mpz_cmp(d,upb)<=0)
    {
      uint64_t i;
      for(i=0;i<cptr;i++)
	{
	  if(mpz_cmp_ui(cs[i],b)<=0)
	    continue;
	  
	  printf("testing r=%lu a=%lu b=%lu c=",ur,a,b);
	  mpz_out_str(NULL,10,cs[i]);
	  printf(" d=");
	  mpz_out_str(NULL,10,d);
	  printf("\n");
	  
	  mpz_mul(tmp1,d,cs[i]);
	  mpz_add_ui(tmp1,tmp1,1);
	  if(mpz_perfect_square_p(tmp1))
	    {
	      	      
	      printf("solution found with r=%lu a=%lu b=%lu c=",ur,a,b);
	      mpz_out_str(NULL,10,cs[i]);
	      printf(" d=");
	      mpz_out_str(NULL,10,d);
	      printf("\n");
	      
	      count++;
	    }
	}
      next_xy(s,t,a,b,ur);
      set_c(s,a,d);
    }
  return(count);
}
*/
mpz_t xn,y_n;

void set_cs_and_ds(uint64_t x, uint64_t y, mpz_t u, mpz_t v, uint64_t D, uint64_t a, uint64_t b1, mpz_t *ztmps)
{
  //printf("in set_cs_and_ds with a=%lu\n",a);
  mpz_set_ui(xn,x);
  mpz_set_ui(y_n,y);
  set_c(xn,a,cs[cptr],lwb,b1); // will set cs[cptr] to 0 if not integral
  //printf("xn=");mpz_out_str(NULL,10,xn);
  //printf(" yn=");mpz_out_str(NULL,10,y_n);
  //printf(" cn=");mpz_out_str(NULL,10,cs[cptr]);printf("\n");
  while(mpz_cmp(cs[cptr],lwb)<=0)
    {
      //printf("in loop\n");
      cptr++;
      mpz_mul(ztmps[0],xn,u);
      mpz_mul_ui(ztmps[2],y_n,D);
      mpz_mul(ztmps[1],ztmps[2],v);
      mpz_mul(ztmps[2],xn,v);
      mpz_mul(ztmps[3],y_n,u);
      mpz_add(xn,ztmps[0],ztmps[1]);
      mpz_add(y_n,ztmps[2],ztmps[3]);
      int res=set_c(xn,a,cs[cptr],lwb,b1);
      if(res==-1)
	{
	  cptr--;
	  break;
	}
      if(res==0)
	cptr--;
    }
  mpz_set(ds[dptr],cs[cptr]);
  while(mpz_cmp(ds[dptr],upb)<=0)
    {
      dptr++;
      mpz_mul(ztmps[0],xn,u);
      mpz_mul_ui(ztmps[2],y_n,D);
      mpz_mul(ztmps[1],ztmps[2],v);
      mpz_mul(ztmps[2],xn,v);
      mpz_mul(ztmps[3],y_n,u);
      mpz_add(xn,ztmps[0],ztmps[1]);
      mpz_add(y_n,ztmps[2],ztmps[3]);
      int res=set_c(xn,a,ds[dptr],upb,b1);
      if(res==-1)
	{
	  dptr--;
	  return;
	}
      if(res==0)
	dptr--;
    }
}

int mycmp(const void *a, const void *b)
{
  return(mpz_cmp(*(mpz_t*)a,*(mpz_t*)b));
}

uint64_t check_ds(mpz_t *ztmps, uint64_t a, uint64_t b)
{
  qsort(cs,cptr,sizeof(mpz_t),mycmp);
  qsort(ds,dptr,sizeof(mpz_t),mycmp);
  uint64_t i,j,count=0;
  for(i=0;i<cptr;i++)
    if(mpz_cmp_ui(cs[i],0)!=0)
      if(mpz_cmp(cs[i],cs[i-1])!=0)
	{
	  mpz_mul(ztmps[0],cs[i],ds[0]);
	  mpz_add_ui(ztmps[0],ztmps[0],1);
	  if(mpz_perfect_square_p(ztmps[0]))
	    {
	      count++;
	      //printf("Solution found %lu %lu ",a,b);
	      //mpz_out_str(NULL,10,cs[i]);printf(" ");
	      //mpz_out_str(NULL,10,ds[0]);printf("\n");
	    }
	  for(j=1;j<dptr;j++)
	    if(mpz_cmp(ds[j],ds[j-1])!=0)
	      {
		mpz_mul(ztmps[0],cs[i],ds[j]);
		mpz_add_ui(ztmps[0],ztmps[0],1);
		if(mpz_perfect_square_p(ztmps[0]))
		  {
		    count++;
		    //printf("Solution found %lu %lu ",a,b);
		    //mpz_out_str(NULL,10,cs[i]);printf(" ");
		    //mpz_out_str(NULL,10,ds[j]);printf("\n");
		  }
	      }
	}
  return(count);
}


uint64_t brute_force(uint64_t a, uint64_t b, mpfi_ptr sqrt_D, mpfi_t *itmps, mpfr_ptr rtmp, mpz_t *ztmps, uint64_t ur)
{
  cptr=0;dptr=0;
  uint64_t g=gcd(b,a),a1=a/g,b1=b/g,D=a1*b1,N=b1*(b1-a1);
  //printf("in brute force with a=%lu b=%lu D=%lu\n",a,b,D);
  //printf("solving X^2-%luy^2==%lu\n",D,N);
  mpz_ui_pow_ui(lwb,b,5);
  mpz_ui_pow_ui(upb,b,8);

  mpfi_set_ui(sqrt_D,D);
  mpfi_sqrt(sqrt_D,sqrt_D);
  abcalls++;
  if(PQm1(pm1res,D,sqrt_D,itmps,rtmp,ztmps[0])==-1)
    {
      //printf("%lu x^2 - %lu y^2 = %lu.\n",b,a,b-a);
      printf("a = %lu b = %lu : X^2 - %lu y^2 = -1 has solution (",a,b,D);
      mpz_out_str(NULL,10,pm1res[0]);printf(",");
      mpz_out_str(NULL,10,pm1res[1]);printf(")\n");

      mpz_mul(ztmps[0],pm1res[0],pm1res[0]);
      mpz_mul(ztmps[2],pm1res[1],pm1res[1]);
      mpz_mul_ui(ztmps[1],ztmps[2],D);
      mpz_mul(pm1res[1],pm1res[0],pm1res[1]);
      mpz_mul_2exp(pm1res[1],pm1res[1],1);
      mpz_add(pm1res[0],ztmps[0],ztmps[1]);
      //printf("X^2-%luy^2=1 has solution (",D);
      //mpz_out_str(NULL,10,pm1res[0]);printf(",");
      //mpz_out_str(NULL,10,pm1res[1]);printf(")\n");

    }


  mpfi_set_ui(itmps[0],N);
  mpz_sub_ui(ztmps[0],pm1res[0],1);
  mpfi_mul_z(itmps[1],itmps[0],ztmps[0]);
  mpfi_div_ui(itmps[0],itmps[1],D<<1);
  mpfi_sqrt(itmps[1],itmps[0]);
  mpfi_floor(ztmps[0],itmps[1],rtmp,ztmps[1]);
  uint64_t y,L2=my_mpz_get_ui(ztmps[0]);
  //printf("L2=%lu\n",L2);
  for(y=0;y<=L2;y++)
    {
      mpz_set_ui(ztmps[0],D);
      mpz_mul_ui(ztmps[1],ztmps[0],y*y);
      mpz_add_ui(ztmps[1],ztmps[1],N);
      if(mpz_perfect_square_p(ztmps[1]))
	{
	  mpz_sqrt(ztmps[0],ztmps[1]);
	  uint64_t x=mpz_get_ui(ztmps[0]); // we know this is <= b
	  if(x%b1==0) // x is an integer
	    {
	      if(y!=1) printf("extra Solution %lu %lu %lu %lu\n",a,b,x,y);
	      set_cs_and_ds(x,y,pm1res[0],pm1res[1],D,a,b1,ztmps);
	      mpz_neg(pm1res[1],pm1res[1]);
	      set_cs_and_ds(x,y,pm1res[0],pm1res[1],D,a,b1,ztmps);
	      mpz_neg(pm1res[1],pm1res[1]);	      
	    }
	}
    }
  int res=check_ds(ztmps,a,b);
  if(res!=2)
    {
      printf("res!=2\n");
    }
  return(res);
}


int main(int argc, char **argv)
{
  if(argc>2)
    {
      printf("Usage:- %s [r factor file].\n",argv[0]);
      exit(0);
    }

  FILE *rfile;
  if(argc==2)
    {
      rfile=fopen(argv[1],"r");
      if(!rfile)
	{
	  printf("Error opening %s for input. Exiting.\n");
	  exit(0);
	}
    }
  else
    rfile=stdin;

  mpfr_set_default_prec(100);
  mpfi_t sqrt_D,itmps[2];
  mpfi_init(sqrt_D);mpfi_init(itmps[0]);mpfi_init(itmps[1]);
  mpfr_t rtmp;
  mpfr_init(rtmp);
  mpz_t ztmps[4];
  mpz_init(ztmps[0]);mpz_init(ztmps[1]);
  
  mpz_init(tmp1);
  mpz_init(tmp2);
  mpz_init(tmp3);
  mpz_init(ab1);
  mpz_init(r);
  mpz_init(s);
  mpz_init(t);
  mpz_init(d);
  mpz_init(upb);
  mpz_init(lwb);
  mpz_init(upbc);
  mpz_init(lwbc);
  mpz_init(rc);
  mpz_init(sc);
  mpz_init(tc);
  mpz_init(c);  
  mpz_init(_P);
  mpz_init(_Q);
  mpz_init(_a);
  mpz_init(_B);
  mpz_init(_B1);
  mpz_init(_B2);
  mpz_init(_G);
  mpz_init(_G1);
  mpz_init(_G2);
  mpz_init(pm1res[0]);
  mpz_init(pm1res[1]);
  mpz_init(xn);
  mpz_init(y_n);

  factor_t f;

  uint64_t i;
  for(i=0;i<C_LIM;i++)
    mpz_init(cs[C_LIM]);
  for(i=0;i<D_LIM;i++)
    mpz_init(ds[D_LIM]);
  for(i=0;i<ZTMPS;i++)
    mpz_init(ztmps[i]);

  //brute_force(21,37,223);exit(0);


  uint64_t r0=0,r1;

  uint64_t ur,uab1,nsols=0;
  while(read_factorisation(rfile,&ur,&f))
    {
      //print_factors(&f,ur*ur-1);
      if(r0==0) r0=ur;
      r1=ur;
      uab1=ur*ur-1;
      uint64_t a=1,b=uab1;
      uint64_t this_pow[MAX_FACS];
      uint64_t i;
      for(i=0;i<f.num_facs;i++)
	this_pow[i]=0;
      while(1==1)
	{
	  if((b<=MAX_B)&&(b>a+2)&&(b<a+a)) nsols+=brute_force(a,b,sqrt_D,itmps,rtmp,ztmps,ur); 
	  this_pow[0]++;
	  uint64_t ptr=0;
	  while(this_pow[ptr]>f.powers[ptr])
	    {
	      this_pow[ptr]=0;
	      a/=f.facs[ptr];
	      ptr++;
	      if(ptr==f.num_facs)
		break;
	      this_pow[ptr]++;
	    }
	  if(ptr==f.num_facs) break;
	  a*=f.primes[ptr];
	  b=uab1/a;
	}
    }

  printf("We found %lu solutions for r in [%lu,%lu].\nWe tried %lu {a,b} pairs\n",nsols,r0,r1,abcalls);
  return(0);
}

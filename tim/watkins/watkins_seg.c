// A segmented version of Watkins' code to confirm L function of
// odd real Dirichlet characters do not have Siegel (exceptional)
// zeros.

#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "math.h"
#include "stdbool.h"

#include "arb.h"


#define INTERESTING_PTR ((int64_t) 0) //290661067L)
#define PRINT_RESULTS (false)
#define VERBOSE (false)

#define PREC (50) // less than double!

int64_t D_max,D_min;

int16_t *counts;
// asqrts[a]=sqrt(a)...
// dlogs[d]=1/2log(d)
// alogs[a]=log(a)
// alphas[d] = sum alpha(a,d)/sqrt(a)
// Hs[d] <= sum H(s;a,b,c)/sqrt(a) for s in [1/2,1)
// S1s[d] = sum (b_3+alpha^3(a,d)/6)/sqrt(a)
// S2s[d] = sum (b_5+alpha^5(a,d)/120)/sqrt(a)

arb_t *asqrts,*dsqrts,*dlogs,*alogs,*alphas,*Hs,*S1s,*S2s;

void S1(arb_t res, arb_t al3)
{
  static bool init=false;
  static arb_t b3,tmp;

  if(!init)
    {
      init=true;
      arb_init(b3);
      arb_init(tmp);
      arb_set_ui(b3,26669630885);
      arb_div_ui(tmp,b3,10,PREC);
      arb_add_error_2exp_si(tmp,-1);
      arb_div_ui(b3,tmp,1000000000,PREC);
      //printf("b3 set to ");arb_printd(b3,20);printf("\n");
    }
  arb_div_ui(tmp,al3,6,PREC);
  arb_add(res,b3,tmp,PREC);
}

void S2(arb_t res, arb_t al5, int64_t d)
{
  static bool init=false;
  static arb_t b5,tmp;

  if(!init)
    {
      init=true;
      arb_init(b5);
      arb_init(tmp);
      arb_set_ui(b5,63999954235);
      arb_div_ui(tmp,b5,10,PREC);
      arb_add_error_2exp_si(tmp,-1); //
      arb_div_ui(b5,tmp,1000000000,PREC);
      //printf("b5 set to ");arb_printd(b5,20);printf("\n");
    }
  arb_div_ui(tmp,al5,120,PREC);
  arb_add(res,b5,tmp,PREC);
}

// compute an upper bound for H(s;a,b,c) for s in [1/2,1)
// see Watkins last display on p 5.
void H(arb_t res, int64_t a, int64_t b, int64_t d, arb_t sqrta)
{
  static bool init=false;
  static arb_t k,two_pi,tmp1,tmp2,tmp3;

  if(!init)
    {
      init=true;
      arb_init(k);
      arb_init(two_pi);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_const_pi(two_pi,PREC);
      arb_mul_2exp_si(two_pi,two_pi,1);
    }

  arb_set_ui(tmp1,b);
  arb_div_ui(tmp2,tmp1,a,PREC);
  arb_cos_pi(tmp1,tmp2,PREC); // cos(pi a/b)
  arb_mul_ui(tmp2,tmp1,50,PREC);
  arb_mul_2exp_si(tmp1,tmp1,1); 
  arb_abs(tmp1,tmp1); // |2cos(pi a/b)|
  arb_add(tmp3,tmp1,tmp2,PREC); // 50 cos - |2 cos|
  
  arb_div_ui(k,dsqrts[d],a+a,PREC); // k=sqrt(d)/2a
  arb_mul(tmp1,k,two_pi,PREC);
  arb_exp(tmp2,tmp1,PREC); // exp(2 pi k)
  arb_mul_ui(tmp1,tmp2,25,PREC);
  arb_sqrt(tmp2,k,PREC);
  arb_mul(k,tmp1,tmp2,PREC); 
  arb_mul(tmp1,k,sqrta,PREC); // 25 sqrt(k) exp(2 pi k) sqrt(a)
  arb_div(res,tmp3,tmp1,PREC);
  
}

// compute alpha(a,d)
void alpha(arb_t res, int64_t a, int64_t d)
{
  static bool init=false;
  static arb_t tmp1,tmp2,log_term;
  if(!init)
    {
      init=true;
      arb_t pi,euler;
      arb_init(pi);
      arb_const_pi(pi,PREC);
      arb_init(euler);
      arb_const_euler(euler,PREC);
      arb_neg(euler,euler);
      arb_init(tmp1);
      arb_exp(tmp1,euler,PREC);
      arb_init(tmp2);
      arb_mul(tmp2,tmp1,pi,PREC);
      arb_mul_2exp_si(tmp2,tmp2,3); // 8 pi exp(-gamma)
      arb_init(log_term);
      arb_log(log_term,tmp2,PREC);
      arb_clear(pi);
      arb_clear(euler);
    }
  arb_add(tmp1,log_term,alogs[a],PREC);
  arb_sub(res,tmp1,dlogs[d],PREC);
}
  

void do_sums(int64_t a, int64_t b, int64_t c, int64_t d, int mult)
{
  static bool init=false;
  static arb_t tmp1,tmp2,al3,al5;
  
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(al3);
      arb_init(al5);
      uint64_t max_a=sqrt((double) D_max/3);
      asqrts=(arb_t *)malloc(sizeof(arb_t)*(max_a+1));
      alogs=(arb_t *)malloc(sizeof(arb_t)*(max_a+1));
      for(uint64_t a=1;a<=max_a;a++)
	{
	  arb_init(asqrts[a]);
	  arb_sqrt_ui(asqrts[a],a,PREC);
	  arb_init(alogs[a]);
	  arb_log_ui(alogs[a],a,PREC);
	}

      arb_t *logs1=(arb_t *)malloc(sizeof(arb_t)*(D_max-D_min+1));
      dlogs=logs1-D_min;
      
      for(uint64_t d=D_min;d<=D_max;d++)
	if(counts[d]>=0)
	  {
	    arb_init(dlogs[d]);
	    arb_log_ui(dlogs[d],d,PREC);
	    arb_mul_2exp_si(dlogs[d],dlogs[d],-1); // contains 1/2 log d
	  }

      logs1=(arb_t *)malloc(sizeof(arb_t)*(D_max-D_min+1));
      alphas=logs1-D_min;
      logs1=(arb_t *)malloc(sizeof(arb_t)*(D_max-D_min+1));
      Hs=logs1-D_min;
      logs1=(arb_t *)malloc(sizeof(arb_t)*(D_max-D_min+1));
      dsqrts=logs1-D_min;
      logs1=(arb_t *)malloc(sizeof(arb_t)*(D_max-D_min+1));
      S1s=logs1-D_min;
      logs1=(arb_t *)malloc(sizeof(arb_t)*(D_max-D_min+1));
      S2s=logs1-D_min;
      for(uint64_t d=D_min;d<=D_max;d++)
	if(counts[d]>=0)
	  {
	    arb_init(alphas[d]);
	    arb_init(Hs[d]);
	    arb_init(dsqrts[d]);
	    arb_sqrt_ui(dsqrts[d],d,PREC);
	    arb_init(S1s[d]);
	    arb_init(S2s[d]);
	  }
    }
  
  alpha(tmp1,a,d);
  arb_div(tmp2,tmp1,asqrts[a],PREC);
  if(mult==2)
    arb_mul_2exp_si(tmp2,tmp2,1);
  arb_add(alphas[d],alphas[d],tmp2,PREC);
  
  arb_mul(tmp2,tmp1,tmp1,PREC); // alpha^2
  arb_mul(al3,tmp2,tmp1,PREC); // alpha^3
  arb_mul(al5,al3,tmp2,PREC); // alpha^5
  
  S1(tmp1,al3);
  arb_div(tmp2,tmp1,asqrts[a],PREC);
  if(mult==2)
    arb_mul_2exp_si(tmp2,tmp2,1);
  arb_add(S1s[d],S1s[d],tmp2,PREC);
  
  S2(tmp1,al5,d);
  arb_div(tmp2,tmp1,asqrts[a],PREC);
  if(mult==2)
    arb_mul_2exp_si(tmp2,tmp2,1);
  arb_add(S2s[d],S2s[d],tmp2,PREC);

  
  H(tmp1,a,b,d,asqrts[a]);
  if(mult==2)
    arb_mul_2exp_si(tmp1,tmp1,1);
  arb_add(Hs[d],Hs[d],tmp1,PREC);
}

void do_checks()
{
  arb_t tmp;
  arb_init(tmp);
  // i) sum alpha(a,d)/sqrt(a) >=0
  // ii) sum(alpha(a,d)-H(s;a,b,c))/sqrt(a) >=0
  for(uint64_t d=D_min;d<=D_max;d++)
    if(counts[d]>=0)
      {
	if(!arb_is_positive(alphas[d]))
	  {
	    printf("Failure in condition 1 with d = -%ld h=%d and sum = ",d,counts[d]);
	    arb_printd(alphas[d],20);
	    printf("\n");
	  }
	arb_sub(tmp,alphas[d],Hs[d],PREC);
	if(!arb_is_positive(tmp))
	  {
	    printf("Failure in condition 2 with d = -%ld h = %d and sum = ",d,counts[d]);
	    arb_printd(tmp,20);
	    printf("\n");
	  }
	if(!arb_is_positive(S1s[d]))
	  {
	    printf("Failure of condition S1 with d = -%ld h = %d and S1 = ",d,counts[d]);
	    arb_printd(S1s[d],20);
	    printf("\n");
	  }
	else
	  {
	    if(!arb_is_positive(S2s[d]))
	      {
		printf("Failure of condition S2 with d = -%ld h = %d and S2 = ",d,counts[d]);
		arb_printd(S2s[d],20);
		printf("\n");
	      }
	  }
      }
  arb_clear(tmp);
    
  return;
}

int main(int argc, char** argv)
{
  printf("Command line:- %s ",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=3)
    {
      printf("Usage:- %s <D_min> <D_max>\n",argv[0]);
      return 0;
    }
  D_min=atol(argv[1]);
  D_max=atol(argv[2]);
  if(D_max<=D_min)
    {
      printf("need max > min.\n");
      return 0;
    }
  if(D_min<4)
    {
      printf("need min > 3.\n");
      return 0;
    }
  if(D_max>10000000000L)
    {
      printf("need max (%ld) <=1e10.\n",D_max);
      return 0;
    }
  
  int16_t *counts1;
  counts1=(int16_t *)malloc(sizeof(int16_t)*(D_max-D_min+1));
  counts=counts1-D_min;

  int64_t a,b,c,A_max,B_max,C_max,ptr,b2,A_min,C_min;

  for(uint64_t i=D_min;i<=D_max;i++)
    {
      uint64_t res16=i&0xF;
      if((res16!=3)&&(res16!=4)&&(res16!=7)&&(res16!=8)&&(res16!=11)&&(res16!=15))
	{
	  counts[i]=-1;
	  continue;
	}
      if((i%9==0)||(i%25==0)||(i%49==0)||(i%121==0)||(i%169)==0)
	{
	  counts[i]=-1;
	  continue;
	}
      counts[i]=0;
    }
  
  // do the b=0 case
  A_max=sqrt((double) D_max);
  A_max>>=1;
  for(a=1;a<=A_max;a++)
    {
      C_min=D_min/a;
      C_min>>=2;
      if(C_min<=a) C_min=a+1;
      ptr=(a*C_min)<<2;
      while(ptr<D_min)
	{
	  C_min++;
	  ptr+=a<<2;
	}	  
      for(c=C_min;;c++)
	{
	  if(ptr>D_max)
	    break;
	  if(INTERESTING_PTR)
	    if(ptr==INTERESTING_PTR)
	      printf("Adding 1, a: %ld b: %ld c:%ld d: %ld\n",a,0l,c,-4*a*c);
	  if(counts[ptr]>=0) {do_sums(a,0,c,ptr,1);counts[ptr]++;}
	  ptr+=a<<2;
	}
    }
  // b=0 case done
  if (VERBOSE) printf("b=0 a!=c done.\n");
  
  B_max=sqrt((double) D_max/(double) 3.0);
  b2=1;
  if(VERBOSE) printf("B_max set to %ld.\n",B_max);
  for(b=1;b<=B_max;b++)
    {
      //b2=b*b;
      // do the a=b case
      if(VERBOSE) if((b&0xff)==1) printf("b = %ld.\n",b);
      a=b;
      C_min=D_min+b2;
      C_min/=(a<<2);
      if(C_min<=a) C_min=a+1;
      ptr=((a*C_min)<<2)-b2;
      while(ptr<D_min)
	{
	  if(VERBOSE) printf("Incrementing C_min (%ld) with a = %ld\n",C_min,a);
	  C_min++;
	  ptr+=a<<2;
	}	  
      for(c=C_min;;c++)
	{
	  if(ptr>D_max)
	    break;
	  if(counts[ptr]>=0) {do_sums(a,b,c,ptr,1);counts[ptr]++;} // a != -b
	  if(INTERESTING_PTR)
	    if(ptr==INTERESTING_PTR)
	      printf("Adding 1, a: %ld b: %ld c: %ld d: %ld\n",a,b,c,b2-4*a*c);
	  ptr+=a<<2;
	}
      // do the a>b case
      A_max=sqrt((double) b2+D_max);
      A_max>>=1;
      for(a=b+1;a<=A_max;a++)
	{
	  C_min=D_min+b2;
	  C_min/=a;
	  C_min>>=2;
	  if(C_min<=a) C_min=a+1;
	  ptr=((a*C_min)<<2)-b2;
	  while(ptr<D_min)
	    {
	      C_min++;
	      ptr+=a<<2;
	    }	  
	  for(c=C_min;;c++)
	    {
	      if(ptr>D_max)
		break;
	      if(counts[ptr]>=0) {do_sums(a,b,c,ptr,2);counts[ptr]+=2;} // +/- b are allowed
	      if(INTERESTING_PTR)
		if(ptr==INTERESTING_PTR)
		  printf("Adding 2, a: %ld b: +/-%ld c: %ld d: %ld\n",a,b,c,b2-4*a*c);
	      ptr+=a<<2;
	    }
	}
      b2+=b+b+1;
    }

  if(VERBOSE) printf("b!=0 c>a done\n");
  
  // do the a=c case;
  // first the b=0 case
  A_min=sqrt((double) D_min);
  A_min>>=2;
  ptr=(A_min*A_min)<<2;
  while(ptr<D_min)
    {
      A_min++;
      ptr+=(A_min+A_min-1)<<2;
    }
  for(a=A_min;;a++)
    {
      if(ptr>D_max)
	break;
      if(counts[ptr]>=0) {do_sums(a,b,a,ptr,1);counts[ptr]++;}
      if(INTERESTING_PTR)
	if(ptr==INTERESTING_PTR)
	  printf("Adding 1, a: %ld b: %ld c: %ld d: %ld\n",a,0l,a,-4*a*a);
      ptr+=(a+a+1)<<2;
    }
  if(VERBOSE) printf("b=0, a=c done.\n");
  b2=1;
  for(b=1;b<=B_max;b++)
    {
      //b2=b*b;
      if(3*b2>D_min)
	if(counts[3*b2]>=0) {do_sums(b,b,b,3*b2,1);counts[3*b2]++;} // a=b=c
      A_min=sqrt((double) D_min+b2);
      A_min>>=1;
      if(A_min<=b) A_min=b+1;
      ptr=((A_min*A_min)<<2)+b2;
      while(ptr<D_min)
	{
	  A_min++;
	  ptr+=(A_min+A_min-1)<<2;
	}
      for(a=A_min;;a++)
	{
	  if(ptr>D_max)
	    break;
	  if(counts[ptr]>=0) {do_sums(a,b,a,ptr,1);counts[ptr]+=1;}
	  if(INTERESTING_PTR)
	    if(ptr==INTERESTING_PTR)
	      printf("Adding 1, a: %ld b: +/-%ld c: %ld d: %ld\n",a,b,a,b2-4*a*a);
	  ptr+=(a+a+1)<<2;
	}
      b2+=b+b+1;
    }
  if(VERBOSE) printf("b!=0 a=c case done.\n");
  if(INTERESTING_PTR)
    printf("%ld %d\n",INTERESTING_PTR,counts[INTERESTING_PTR]);

  if(VERBOSE)
    for(int64_t d=D_min;d<=D_max;d++)
      if(counts[d]>=0)
	{printf("d: %ld h: %d alpha_sum: ",d,counts[d]);arb_printd(alphas[d],20);printf(" Hs: ");arb_printd(Hs[d],20);printf("\n");}


  do_checks();
  
  return 0;
}

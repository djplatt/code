// A segmented version of Watkins' code to confirm L function of
// odd real Dirichlet characters do not have Siegel (exceptional)
// zeros.
// uses the int_double interval arithmetic package

#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "math.h"
#include "stdbool.h"
#include "primesieve.h"

#include "int_double12.0.h"

#define VERBOSE (false)

int64_t D_max,D_min;

// counts[d]=class number of Q[-d]
// asqrts[a]=sqrt(a)
// dlogs[d]=1/2log(d)
// alogs[a]=log(a)
// S0s[d] = sum_Q alpha(a,d)/sqrt(a)
// Hs[d] <= sum_Q H(s;a,b,c)/sqrt(a) for s in [1/2,1)
// S1s[d] = sum_Q (b_3+alpha^3(a,d)/6)/sqrt(a)
// S2s[d] = sum_Q (b_5+alpha^5(a,d)/120)/sqrt(a)

int64_t *counts;
int_double *asqrts,*dsqrts,*dlogs,*alogs,*S0s,*Hs,*S1s,*S2s;
int_double *tmp1,*tmp2,*tmp3,*tmp4,*tmp5,*tmp6;

int_double S1(int_double al3)
{
  static bool init=false;
  static int_double b3;

  if(!init)
    {
      init=true;
      b3=int_double(2666963088,2666963089);
      b3/=1000000000;
    }
  return b3+al3/6;
}

int_double S2(int_double al5, int64_t d)
{
  static bool init=false;
  static int_double b5;

  if(!init)
    {
      init=true;
      b5=int_double(6399995423,6399995424);
      b5/=1000000000;
    }
  
  int_double res=b5+al5/120;
  return res;
}

// in here we have x = b/a so |x|<=1
int_double my_cospi(int_double x)
{
  int_double pi_x2=sqr(d_pi*x);
  int_double res=1-pi_x2/2.0;
  int_double pi_xn=sqr(pi_x2); // (pi x)^4
  res+=pi_xn/24.0;
  pi_xn*=pi_x2; // (pi x)^6
  res-=pi_xn/720.0;
  pi_xn*=pi_x2; // (pi x)^8
  res+=pi_xn/40320.0;
  /*
  int_double cos_err=pi_xn*pi_x2/3628800.0; // [e,e]
  cos_err.left=0.0; // [0,e]
  return res-cos_err;
  */
  pi_xn*=pi_x2; // (pi x)^10
  res-=pi_xn/3628800.0;
  pi_xn*=pi_x2; // (pi x)^12
  res+=pi_xn/479001600.0;
  pi_xn*=pi_x2; // (pi x)^14
  res-=pi_xn/87178291200.0;
  return res+int_double(0.0,4.3031e-6); // [0,pi^16/16!]
}
  

// compute an upper bound for H(s;a,b,c)/sqrt(a) for s in [1/2,1)
// see Watkins last display on p 5.
/*int_double H(int64_t d)
{
  static bool init=false;
  static int_double c;
  if(!init)
    {
      init=true;
      c=127471;
      c/=10000000;
    }
  return c*dqrts[d];
}
*/

int_double H(int64_t a, int64_t b, int64_t d, int_double sqrta)
{
  int_double k=sqrt(int_double(d))/(2*a);
  int_double cos_bit=my_cospi(int_double(b)/a);
  int_double exp_bit=exp(d_two_pi*k);

  int_double res=(50*cos_bit+abs(cos_bit))/(25*sqrt(k)*exp_bit*sqrta);
  
  //if(d==9999999043) print_int_double_str("H returning ",res);
  return res;
}


// compute alpha(a,d)
int_double alpha(int64_t a, int64_t d)
{
  static bool init=false;
  static int_double log_term;
  if(!init)
    {
      init=true;
      int_double d_gamma=int_double(577215664901,577215664902);
      d_gamma/=1000000000000L; // Euler's gamma = 0.577...
      log_term=log(8*d_pi*exp(-d_gamma));
    }
  return log_term+alogs[a]-dlogs[d];
}

// we have found -d=b^2-4ac meeting the criteria
// if mult = 2 then -b works as well
void do_sums(int64_t a, int64_t b, int64_t c, int64_t d, int mult)
{
  static bool init=false;
  
  if(!init) // compute tables once and for all
    {
      init=true;
      uint64_t max_a=sqrt((double) D_max/3);
      asqrts=(int_double *)malloc(sizeof(int_double)*(max_a+1));
      alogs=(int_double *)malloc(sizeof(int_double)*(max_a+1));
      for(uint64_t a=1;a<=max_a;a++)
	{
	  int_double aa=a;
	  asqrts[a]=sqrt(aa);
	  alogs[a]=log(aa);
	}

      tmp1=(int_double *)malloc(sizeof(int_double)*(D_max-D_min+1));
      dlogs=tmp1-D_min;
      tmp2=(int_double *)malloc(sizeof(int_double)*(D_max-D_min+1));
      dsqrts=tmp2-D_min;
      tmp3=(int_double *)malloc(sizeof(int_double)*(D_max-D_min+1));
      S0s=tmp3-D_min;
      tmp4=(int_double *)malloc(sizeof(int_double)*(D_max-D_min+1));
      Hs=tmp4-D_min;
      tmp5=(int_double *)malloc(sizeof(int_double)*(D_max-D_min+1));
      S1s=tmp5-D_min;
      tmp6=(int_double *)malloc(sizeof(int_double)*(D_max-D_min+1));
      S2s=tmp6-D_min;
      for(uint64_t d=D_min;d<=D_max;d++)
	if(counts[d]>=0)
	  {
	    int_double dd=d;
	    dlogs[d]=log(dd)/2.0;
	    dsqrts[d]=sqrt(dd);
	  }
    }

  if(VERBOSE)
    {
      if(mult==2)
	printf("doing a: %ld +/-b: %ld c: %ld d: %ld\n",a,b,c,d);
      else
	printf("doing a: %ld b: %ld c: %ld d: %ld\n",a,b,c,d);
    }
  
  int_double alph=alpha(a,d);
  if (VERBOSE) {printf("alpha returned ");print_int_double(alph);printf("\n");}
  if(mult==2)
    S0s[d]+=(alph+alph)/asqrts[a];
  else
    S0s[d]+=alph/asqrts[a];
  
  int_double al3=alph*alph*alph;
  int_double al5=al3*alph*alph;
  
  int_double this_S1=S1(al3);
  if(VERBOSE) {printf("S1 returned ");print_int_double(this_S1);printf("\n");}
  if(mult==2)
    S1s[d]+=(this_S1+this_S1)/asqrts[a];
  else
    S1s[d]+=this_S1/asqrts[a];

  int_double this_S2=S2(al5,d);
  if(VERBOSE) {printf("S2 returned ");print_int_double(this_S2);printf("\n");}
  if(mult==2)
    S2s[d]+=(this_S2+this_S2)/asqrts[a];
  else
    S2s[d]+=this_S2/asqrts[a];

  int_double this_H=H(a,b,d,asqrts[a]);
  if(VERBOSE) {printf("H returned ");print_int_double(this_H);printf("\n");}

  if(mult==2)
    Hs[d]+=this_H+this_H;
  else
    Hs[d]+=this_H;
}


void do_checks()
{
  // i) sum alpha(a,d)/sqrt(a) >=0
  // ii) sum (alpha(a,d)-H(s;a,b,c))/sqrt(a) >=0
  // iii) sum S1 >=0
  for(uint64_t d=D_min;d<=D_max;d++)
    if(counts[d]>=0)
      {
	if(S0s[d].left<0.0)
	  {
	    printf("Failure in condition 1 with d = -%ld h = %ld and sum = ",d,counts[d]);
	    print_int_double(S0s[d]);
	    printf("\n");
	  }
	int_double tmp=S0s[d]-Hs[d];
	if(tmp.left<0.0)
	  {
	    printf("Failure in condition 2 with d = -%ld h = %ld and sum = ",d,counts[d]);
	    print_int_double(tmp);
	    printf("\n");
	  }
	if(S1s[d].left<0.0)
	  {
	    printf("Failure of condition S1 with d = -%ld h = %ld and S1 = ",d,counts[d]);
	    print_int_double(S1s[d]);
	    printf("\n");
	  }
	if(S2s[d].left<0.0)
	  {
	    printf("Failure of condition S2 with d = -%ld h = %ld and S2 = ",d,counts[d]);
	    print_int_double(S2s[d]);
	    printf("\n");
	  }
      }
    
  return;
}

void do_hs()
{
  int64_t max_h=0;
    for(uint64_t d=D_min+1;d<=D_max;d++)
      if(counts[d]>max_h)
	max_h=counts[d];
    int64_t *hs=(int64_t *) malloc(sizeof(int64_t)*(max_h+1));
    int64_t *last_hs=(int64_t *) malloc(sizeof(int64_t)*(max_h+1));
    for(int64_t h=0;h<=max_h;h++)
      hs[h]=0;
    for(uint64_t d=D_min+1;d<=D_max;d++)
      {
	int64_t h=counts[d];
	if(h>=0)
	  {
	    last_hs[h]=d;
	    if(hs[h]==0)
	      printf("First h = %ld is at d = -%lu.\n",h,d);
	    hs[h]++;
	  }
      }
    for(int64_t h=0;h<=max_h;h++)
      if(hs[h]>0)
	printf("There are %ld discriminants with class number %ld.\n",hs[h],h);
    for(int64_t h=1;h<=max_h;h++)
      if(last_hs[h]>0)
	printf("Largest d with class number %ld was %ld.\n",h,last_hs[h]);    
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
  if(D_min<=4)
    {
      printf("need min > 4.\n");
      return 0;
    }
  if(D_max>10000000000L)
    {
      printf("need max (%ld) <=1e10.\n",D_max);
      return 0;
    }

  _fpu_rndd();
  
  int64_t *counts1;
  counts1=(int64_t *)malloc(sizeof(int64_t)*(D_max-D_min+1));
  counts=counts1-D_min;

  int64_t a,b,c,A_max,B_max,C_max,ptr,b2,A_min,C_min;

  // remove (most) non-FD
  for(uint64_t i=D_min;i<=D_max;i++)
    {
      counts[i]=0;
      uint64_t res16=i&0xF;
      if((res16!=3)&&(res16!=4)&&(res16!=7)&&(res16!=8)&&(res16!=11)&&(res16!=15))
	counts[i]=-1;
    }

  // cross out multiples of odd prime squares
  primesieve_iterator it;
  primesieve_init(&it);
  uint64_t p=primesieve_next_prime(&it);
  while(true)
    {
      p=primesieve_next_prime(&it); // start at 3
      uint64_t p2=p*p;
      if(p2>D_max)
	break;
      ptr=(D_min/p2)*p2;
      while(ptr<D_min)
	ptr+=p2;
      while(ptr<=D_max)
	{
	  counts[ptr]=-1;
	  ptr+=p2;
	}
    }
  primesieve_free_iterator(&it);
  
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
	  if(counts[ptr]>=0)
	    {
	      do_sums(a,0,c,ptr,1);
	      counts[ptr]++;
	    }
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
	  if(counts[ptr]>=0)
	    {
	      do_sums(a,b,c,ptr,1);
	      counts[ptr]++;
	    } // a != -b
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
	      if(counts[ptr]>=0)
		{
		  do_sums(a,b,c,ptr,2);
		  counts[ptr]+=2;
		} // +/- b are allowed
	      ptr+=a<<2;
	    }
	}
      b2+=b+b+1;
    }

  if(VERBOSE) printf("b!=0 c>a done\n");
  
  // do the a=c case;
  b2=0;
  for(b=0;b<=B_max;b++)
    {
      A_min=sqrt((double) D_min+b2);
      A_min>>=1;
      if(A_min<b)
	A_min=b;
      ptr=((A_min*A_min)<<2)-b2;
      while(ptr<D_min)
	{
	  A_min++;
	  ptr+=(A_min+A_min-1)<<2;
	}
      for(a=A_min;;a++)
	{
	  if(ptr>D_max)
	    break;
	  if(counts[ptr]>=0)
	    {
	      do_sums(a,b,a,ptr,1);
	      counts[ptr]++;
	    }
	  ptr+=(a+a+1)<<2;
	}
      b2+=b+b+1;
    }

  if(VERBOSE) printf("a=c case done.\n");

  do_checks();
  free(asqrts);
  free(alogs);
  free(tmp1);
  free(tmp2);
  free(tmp3);
  free(tmp4);
  free(tmp5);
  free(tmp6);

  do_hs();
  
  return 0;
}

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define MAX_FACS (16) // product of 1st 16 primes > 2^64
// structure to hold factor information
typedef struct
{
        unsigned long int pr;   // = 0 if no primitive root
        unsigned long int phi;  // Euler phi(q)
        unsigned long int num_facs;  // number of prime power factors
        unsigned long int facs[MAX_FACS];  // the factors p^n
        unsigned long int primes[MAX_FACS]; // the prime factors p
} factor;

long unsigned int gcd (long unsigned int a, long unsigned int b)
/* Euclid algorithm gcd */
{
	long unsigned int c;
	while(a!=0)
	{
		c=a;
		a=b%a;
		b=c;
	};
	return(b);
};

long unsigned int co_prime(long unsigned int a, long unsigned int b)
{
	return(gcd(a,b)==1);
};


long unsigned int pow_mod(long unsigned int a, long unsigned int b, long unsigned int m)
{
  long unsigned int a_pw=a,pw=b,res=1;
  while(true)
    {
      //printf("pw=%ld a_pw=%ld res=%ld\n",pw,a_pw,res);
      if(pw&1)
	res=(res*a_pw)%m;
      pw>>=1;
      if(pw==0)
	return(res);
      a_pw=(a_pw*a_pw)%m;
    }
}

// find a primitive root for i
// assumes there is one
unsigned long int pr(unsigned long int i, factor *factors)
{
  unsigned long int phi=factors[i].phi;
  for(unsigned long int p=2;p<i;p++)
    {
      if(gcd(p,i)!=1)
	continue;
      bool good=true;
      for(unsigned long int j=0;j<factors[phi].num_facs;j++)
	{
	  if(pow_mod(p,phi/factors[phi].primes[j],i)!=1)
	    continue;
	  good=false;
	  break;
	}
      if(good)
	return(p);
    }
}

bool make_factors(factor *factors, unsigned long int q_end)
{
  if(!factors)
    return false;
  unsigned long int *primes,max_p=floor(sqrt(q_end)),f,i,n,j,p2;
  primes=(unsigned long int *) malloc(sizeof(unsigned long int)*(q_end+1));
  for(i=2;i<=q_end;i++)
    primes[i]=0;
  for(i=4;i<=q_end;i+=2)
    primes[i]=2;
  for(i=3;i<=max_p;i+=2)
    if(primes[i]==0)
      for(j=i*i;j<=q_end;j+=i)
	if(primes[j]==0)
	  primes[j]=i;
  printf("Prime sieve completed.\n");
  // now each entry primes[i] is 0 if i prime, else = smallest prime factor
  factors[3].num_facs=1;
  factors[3].primes[0]=3;
  factors[3].facs[0]=3;
  factors[4].num_facs=1;
  factors[4].primes[0]=2;
  factors[4].facs[0]=4;

  for(f=5;f<=q_end;f++)
    if(primes[f]==0) // a prime
      {
	factors[f].num_facs=1;
	factors[f].primes[0]=f;
	factors[f].facs[0]=f;
      }
    else
      {
	factors[f].primes[0]=primes[f];
	n=f/primes[f];
	if(factors[n].primes[0]==primes[f])
	  {
	    factors[f].num_facs=factors[n].num_facs;
	    factors[f].facs[0]=primes[f]*factors[n].facs[0];
	    for(i=1;i<factors[n].num_facs;i++)
	      {
		factors[f].primes[i]=factors[n].primes[i];
		factors[f].facs[i]=factors[n].facs[i];
	      }
	  }
	else
	  {
	    factors[f].num_facs=factors[n].num_facs+1;
	    factors[f].facs[0]=primes[f];
	    factors[f].primes[0]=primes[f];
	    for(i=1;i<factors[f].num_facs;i++)
	      {
		factors[f].primes[i]=factors[n].primes[i-1];
		factors[f].facs[i]=factors[n].facs[i-1];
	      }
	  }
      }
  free(primes);
  printf("Factors computed.\n");
  // now calculate phi(f)
  factors[1].phi=1;
  for(i=2;i<=q_end;i++)
    {
      if(factors[i].num_facs==1) // prime or prime power
	{
	  factors[i].phi=factors[i].facs[0]-factors[i].facs[0]/factors[i].primes[0]; // p^n - p^(n-1)
	  continue;
	}
      // its a product of prime powers, all of which will be known
      factors[i].phi=1;
      for(j=0;j<factors[i].num_facs;j++)
	factors[i].phi*=factors[factors[i].facs[j]].phi;
    }
  printf("phi computed.\n");

  //now do the prim roots
  //q_end=50;printf("q_end=%ld\n",q_end);
  for(i=3;i<=q_end;i++)
    factors[i].pr=0;
  factors[4].pr=3;
  for(i=3;i<=q_end;i++)
    {
      if(factors[i].pr!=0) // already found one
	{
	  //printf("%ld already has a primitive root, skipping.\n",i);
	  continue;
	}
      if(factors[i].primes[0]==i) // its an odd prime
	{
	  //printf("Doing i=%ld\n",i);
	  factors[i].pr=pr(i,factors); // the base case
	  //printf("pr(%ld) set to %ld\n",i,factors[i].pr);
	  p2=i*i;
	  if(p2<=q_end)
	    {
	      if(pow_mod(factors[i].pr,i-1,p2)==1)
		for(j=p2;j<=q_end;j*=i)
		  {
		    factors[j].pr=factors[i].pr+i;
		    //printf("pr(%ld) set to %ld\n",j,factors[j].pr);
		  }
	      else
		for(j=p2;j<=q_end;j*=i)
		  {
		    factors[j].pr=factors[i].pr;
		    //printf("pr(%ld) set to %ld\n",j,factors[j].pr);
		  }
	    }
	  p2=i<<1;
	  if(p2<=q_end)
	    {
	      if((factors[i].pr&1)==1)
		for(j=p2;j<=q_end;j*=i)
		  factors[j].pr=factors[i].pr;
	      else
		for(j=p2;j<=q_end;j*=i)
		  factors[j].pr=factors[i].pr+(j>>1);
	    }
	  continue;
	}
      factors[i].pr=0;
    }
  printf("pr's computed.\n");
  return(true);
}


// on exit q_n[(pr^n) mod q]=n
void fill_q_ns (unsigned long int *q_n,
		unsigned long int pr, unsigned long int phi_q,
		unsigned long int q)
{
	unsigned long int i,pr_n;
	q_n[1]=0;
	q_n[pr]=1;
	pr_n=pr;
	for(i=2;i<phi_q;i++)
	{
		pr_n=(pr_n*pr)%q;
		q_n[pr_n]=i;
	}
}


void fill_q_ns_2s(unsigned long int *q_n, unsigned long int q)
{
	unsigned long int i,pr;
	pr=1;
	for(i=0;i<(q>>2);i++)
	{
		q_n[pr]=i;
		q_n[q-pr]=i+(q>>2);
		pr=(pr*5)%q;
	}
}

// on exit, offset[i] is where the i'th fraction a/q goes in the fft vector
void make_offsets(unsigned long int q, unsigned long int *q_n,
		  unsigned long int *offset_vec, factor *factors,
		  unsigned long int *a_n)
{
  unsigned long int i,j,fac,fac1,offset,offset1,ptr;

	ptr=0;   // ptr into offset vec
	fac=factors[q].facs[0];
	for(i=0;i<factors[q].phi;i++)
	{
		offset=q_n[a_n[i]%fac];
		offset1=fac;
		for(j=1;j<factors[q].num_facs;j++)
		{
			fac1=factors[q].facs[j];
			offset*=factors[fac1].phi;
			offset+=q_n[a_n[i]%fac1+offset1];
			offset1+=fac1;
		}
		offset_vec[ptr++]=offset;
	}
}


int main(int argc, char** argv)
{
  if(argc!=2)
    return 0;
  unsigned long int q = atoi(argv[1]);
  if((q%4)==2) // no primitive characters
    return 0;
  
  factor *factors=(factor *)malloc(sizeof(factor)*(q+1));
  if(!make_factors(factors,q))
    {
      printf("Failed to construct factors database.\n");
      return 0;
    }

  unsigned long int *a_n=(unsigned long int *)malloc(sizeof(unsigned long int)*factors[q].phi);

  unsigned long int i,j;
  j=0;
  for(i=1;i<q;i++)
    if(co_prime(i,q))
      a_n[j++]=i;

  
  unsigned long int *q_n=(unsigned long int *)malloc(sizeof(unsigned long int)*q);
  unsigned long int *offset_vec=(unsigned long int *)malloc(sizeof(unsigned long int)*q);

  if(factors[q].pr!=0) // q has a primitive root, so nice easy 1-d fft
    {
      fill_q_ns(q_n,factors[q].pr,factors[q].phi,q);

      for(i=1;i<q;i++)
	if(co_prime(i,q))
	  printf("Put zeta(1/2,%lu/%lu) in location %lu.\n",i,q,q_n[i]);
      return 0;
    }      

      //printf("%ld does not have a generator.\n",q);
  long unsigned int no_dims=factors[q].num_facs;
  long unsigned int dims[MAX_FACS+1]; // or malloc no_dims+1
  long unsigned int fac=factors[q].facs[0];        // the first p^n
  bool power_2=(factors[fac].pr==0);  // no conductor => p=2, n>=3
  if(power_2)
    {
      //printf("%ld contains a power of 2^n, n>=3\n");
      no_dims++;
      fill_q_ns_2s(q_n,fac);    // use the {-1,1}X{5} trick
      for(i=1;i<factors[q].num_facs;i++) // move all the dimensions up one
	dims[i+1]=factors[factors[q].facs[i]].phi;
      dims[1]=factors[q].facs[0]/4;
      dims[0]=2;                         // slip in a two
    }
  else
    {
      //printf("Filling q_ns for factor %ld.\n",fac);
      fill_q_ns(q_n,factors[fac].pr,factors[fac].phi,fac); // use the generator
      //printf("q_ns filled.\n");
      for(i=0;i<factors[q].num_facs;i++)
	dims[i]=factors[factors[q].facs[i]].phi;
    }
  long unsigned int offset=fac;
  for(i=1;i<factors[q].num_facs;i++)  // do the rest on the factors, all will have generators
    {
      fac=factors[q].facs[i];      
      //printf("Filling q_ns for factor %ld.\n",fac);
      fill_q_ns(&q_n[offset],factors[fac].pr,factors[fac].phi,fac);  // use the generator
      //printf("q_ns filled.\n");
      offset+=fac;
    }
  //printf("making offsets.\n");
  make_offsets(q,q_n,offset_vec,factors,a_n);    // reverse q_n so we know how to populate fft vector
  for(i=0;i<factors[q].phi;i++)
    printf("Put zeta(1/2,%lu/%lu) in %lu.\n",a_n[i],q,offset_vec[i]);

  
  return 0;
}

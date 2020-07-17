#ifndef MAKE_FACTORS
#define MAKE_FACTORS

#ifndef TRUE
#define TRUE (1==1)
#define FALSE (1==0)
#endif

inline long int phi(long int p, long int p_n)
{
  return(p_n-p_n/p);
}

inline long int mod_pow(long int a, long int b, long int m) // a^b mod m
{
  //printf("%ld^%ld mod %ld returning ",a,b,m);

  long int res=1,pow=a;

  while(b>0)
    {
      if(b&1)
	res=(res*pow)%m;
      pow=(pow*pow)%m;
      b>>=1;
    }
  //printf("%ld\n",res);
  return(res);
}


long int pr(long int n, long int phi_n, factor *factors)
{
  long int pr,fac;
  int good;
  //printf("in pr with n=%ld phi(n)=%ld\n",n,phi_n);
  for(pr=2;;pr++) // there better be one
    {
      /*
      if(n==961)
	printf("trying pr=%ld\n",pr);
      */
      if(gcd(n,pr)!=1)
	continue;
      good=TRUE;
      for(fac=0;fac<factors[phi_n].num_facs;fac++)
	if(mod_pow(pr,phi_n/factors[phi_n].primes[fac],n)==1)
	  {
	    good=FALSE;
	    break;
	  }
      if(good)
	{
	  //printf("pr returning %ld\n",pr);
	  return(pr);
	}

    }
}
	

int make_factors(factor *factors, long int q_end)
{
  long int *primes,max_p=floor(sqrt(q_end)),f,i,n,j,p2;
  primes=(long int *) malloc(sizeof(long int)*(q_end+1));
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
  for(i=3;i<=q_end;i++)
    {
      factors[i].phi=1;
      for(j=0;j<factors[i].num_facs;j++)
	factors[i].phi*=phi(factors[i].primes[j],factors[i].facs[j]);
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
	  factors[i].pr=pr(i,factors[i].phi,factors); // the base case
	  //printf("pr(%ld) set to %ld\n",i,factors[i].pr);
	  p2=i*i;
	  if(p2<=q_end)
	    {
	      if(mod_pow(factors[i].pr,i-1,p2)==1)
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
  /*
  for(int i=3;i<50;i++)
    {
      printf("%ld has %ld factors, phi(%ld)=%ld, pr(%ld)=%ld, factors are",i,factors[i].num_facs,i,factors[i].phi,i,factors[i].pr);
      for(int j=0;j<factors[i].num_facs;j++)
	printf(" %ld %ld",factors[i].primes[j],factors[i].facs[j]);
      printf("\n");
    }
  exit(0);
  */
  return(TRUE);
}

#endif

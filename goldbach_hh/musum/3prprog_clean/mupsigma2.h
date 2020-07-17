typedef struct primfact {
  long int p;
  short int r;
  primfact *next;
} primfact;

#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

void fillisprime(short int *isprime, long N)
/* sets isprime[n]=0 if n is composite, isprime[n]=1 i n is prime
   for 2<=n<N */
{
  long i,j;

  for(i=2; i<N; i++)
    isprime[i]=1;
  for(i=2; i*i<N; i++)
    if(isprime[i])
      for(j=i; j*i<N; j++)
        isprime[j*i]=0;
}

void fillmublock(short int *mun, short int *isprime,
		 long int n, long int m)
/* fills mun[0], mun[1],... with mu(n), mu(n+1), ..., mu(n+m-1) */
/* assumes isprime is filled and valid up to and including sqrt(n+m) */
/* convention: mu(n)=0 */
{
  long int i,j;
  long int *tabla;
  long int maxp;
  long int p, red;

  tabla = (long int *) calloc(m,sizeof(long int));
  for(i=0; i<m; i++)
    tabla[i] = 1;
  
  for(i=0; i<m; i++) {
    if((i+n)%2)
      mun[i]=1;
    else if((i+n)%4) {
      mun[i]=-1;
      tabla[i] = 2;
    }
    else mun[i]=0;
  }
  
  for(p=3; p*p<=n+m; p++) 
    if(isprime[p]) {
      red = (n%p);
      for(j= (red? p-red : 0); j<m; j+=p) {
	if(((j+n)/p)%p) {
	  mun[j] = - mun[j];
	  tabla[j] *= p;
	}
	else
	  mun[j] = 0;
      }
    }
  
  for(i=0; i<m; i++)
    if(mun[i]!=0 && tabla[i]!=n+i)
      mun[i] = - mun[i];

  free(tabla);
}


void fillcopblock(short int *cun, short int *isprime,
		  long int n, long int m, long int r)
/* fills cun[0],...,cun[m-1] as follows:
    cun[i] is 1 if n+i and r are coprime, 
    cun[i] is 0 otherwise */
/* assumes isprime is filled and valid up to and including sqrt(n+m) */
/* assumes r>0 */
/* convention: mu(n)=0 */

/*very stupid algorithm; using a factorization of r would be better */
{
  long int p, i, j, nr, red;
  
  for(j=0; j<=m-1; j++)
    cun[j]=1;

  for(p=2; p*p<=r; p++) 
    if(isprime[p] && !(r%p)) {
      red = (n%p);
      for(j= (red? p-red : 0); j<m; j+=p) 
	cun[j] = 0;
      do {r/=p;} while(!(r%p));
    }

  if(r>1) {
    nr = n%r;
    for(j= (nr ? r- nr : 0); j<m; j+=r) 
      cun[j] = 0;
  }
}

void fillsigmablock(long int *sun, short int *isprime,
		    long int n, long int m)
/* fills sun[0], sun[1],... 
 with sigma(n), sigma(n+1), ..., sigma(n+m-1) */
/* assumes isprime is filled and valid up to and including sqrt(n+m) */
/* gives correct results only for n such that mu(n)=0 */
{
  long int i,j;
  long int *tabla;
  long int p, red;

  tabla = (long int *) calloc(m,sizeof(long int));
 
  for(i=0; i<m; i++) {
    if((i+n)%2) {
      sun[i]=1; tabla[i] = 1;
    }  else {
      sun[i]=3; tabla[i] = 2;
    }
  }

  for(p=3; p*p<=n+m; p+=2) 
    if(isprime[p]) {
      red = (n%p);
      for(j= (red? p-red : 0); j<m; j+=p) {
	sun[j] *= p+1;
	tabla[j] *= p;
      }
    }

  for(i=0; i<m; i++)
    if(tabla[i]!=n+i)
      sun[i] *= ((n+i)/tabla[i]) + 1;

  free(tabla);
}

void fillphiblock(long int *phun, short int *isprime,
		    long int n, long int m)
/* fills phun[0], phun[1],... 
 with phi(n), phi(n+1), ..., phi(n+m-1) */
/* assumes isprime is filled and valid up to and including sqrt(n+m) */
/* gives correct results only for n such that mu(n)=0 */
{
  long int i,j;
  long int *tabla;
  long int p, red;

  tabla = (long int *) calloc(m,sizeof(long int));
 
  for(i=0; i<m; i++) {
    if((i+n)%2) {
      phun[i]=1; tabla[i] = 1;
    }  else {
      phun[i]=1; tabla[i] = 2;
    }
  }

  for(p=3; p*p<=n+m; p+=2) 
    if(isprime[p]) {
      red = (n%p);
      for(j= (red? p-red : 0); j<m; j+=p) {
	phun[j] *= p-1;
	tabla[j] *= p;
      }
    }

  for(i=0; i<m; i++)
    if(tabla[i]!=n+i)
      phun[i] *= ((n+i)/tabla[i]) - 1;

  free(tabla);
}


void fillfactblock(primfact **fun, short int *isprime, long int n, long int m) 
/* fills fun[0], fun[1],... 
 with the list of prime factors of sigma(n), sigma(n+1), ..., sigma(n+m-1) */
/* assumes isprime is filled and valid up to and including sqrt(n+m) */
/* convention: 0 and 1 factorize into nothing */
{
  long i,j,j0;
  primfact *ncell;
  long int *tabla;
  long int p;

  tabla = (long int *) calloc(m,sizeof(long int));
  for(i=0; i<m; i++)
    tabla[i] = n+i;
  memset(fun, 0, sizeof(primfact *)*m);
  
  for(p=2; p*p<=n+m; p++) 
    if(isprime[p]) {
          if(n==0) j0=p; else j0 = ((n+p-1)/p)*p-n;
      for(j=j0; j<m; j+=p) {
	ncell = (primfact *) malloc(sizeof(primfact));
	ncell->p = p; ncell->next = fun[j]; ncell->r=0;
	do {
	  tabla[j] /= p; ncell->r++;
	} while(!(tabla[j]%p));
	fun[j]=ncell;
      }
    }

  for(i=0; i<m; i++)
    if(tabla[i]>1) {
      ncell = (primfact *) malloc(sizeof(primfact));
      ncell->p = tabla[i]; ncell->next = fun[i]; ncell->r=1;
      fun[i]=ncell;
    }

  free(tabla);
}


long int *subdivlist(primfact *fl, long int *oul, long int prod);

long int *divlist(primfact *fl)
/* given the list of prime divisors of a number, return a list containing
   all divisors, terminated by 0*/
{
  primfact *pt;
  long int *oul;
  long len;

  for(len=1, pt=fl; pt; pt = pt->next)
    len *= (pt->r+1);

  oul = (long int *)  calloc(len+1,sizeof(long int));
  *(subdivlist(fl,oul,1)) = 0;

  return oul;
}

long int *subdivlist(primfact *fl, long int *oul, long int prod)
{
  short int i;
  
  if(!fl) {
    *oul++ = prod;
    return oul;
  }

  for(i=0; i<=fl->r; i++) {
    oul=subdivlist(fl->next,oul,prod);
    prod *= fl->p;
  }
  
  return oul;
}

long int *subdivsiglist(primfact *fl, long int **divo, long int **sigo,
			long int prod, long int sprod);

long int *divsiglist(primfact *fl, long int **divo, long int **sigo)
/* given the list of prime divisors of a number, return a list containing
   all divisors, terminated by 0*/
{
  primfact *pt;
  long len;
  long int *doul, *soul;
  
  for(len=1, pt=fl; pt; pt = pt->next)
    len *= (pt->r+1);

  *divo = doul = (long int *)  calloc(len+1,sizeof(long int));
  *sigo = soul = (long int *)  calloc(len+1,sizeof(long int));

  subdivsiglist(fl,divo,sigo,1,1);
  **divo = 0; **sigo = 0;

  *divo = doul; *sigo = soul;
}

long int *subdivsiglist(primfact *fl, long int **divo, long int **sigo,
			long int prod, long int sprod)
{
  short int i;
  
  if(!fl) {
    *(*divo)++ = prod;
    *(*sigo)++ = sprod;
  } else
    for(i=0; i<=fl->r; i++) {
      subdivsiglist(fl->next,divo,sigo,prod,sprod);
      prod *= fl->p;
      if(!i)
	sprod *= fl->p+1;
      else
	sprod *= fl->p;
    }
}


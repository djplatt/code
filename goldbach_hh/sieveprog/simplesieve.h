short *smallprime=NULL;
unsigned long smallM;

void SimpleSiev(short *p, unsigned long N)
/* assumes an array of length >=N+1 is allocated at p */
/* ensures: for 1<=n<=N, p[n]=1 if n is prime, p[n}=0 otherwise*/
{
  unsigned long n,m;

  if(N<1) return;
  p[1]=0;
  if(N<2) return;
  p[2]=1;
  if(N<3) return;
  
  for(n=3; n<=N-1;) {
    p[n++]=1; p[n++]=0;
  }
  if(n==N)
    p[n]=1;

  m=3; n=m*m;
  for(; n<=N;) {
    for(; n<=N ; n+=2*m)
      p[n]=0;
    m+=2; n = m*m;
  }
}

void AllocSimpleSiev(unsigned long M)
{
  if(!smallprime || smallM<M) {
    if(smallprime) free(smallprime);
    smallprime = (short *) calloc(M+1,sizeof(short));
    smallM= M;
    SimpleSiev(smallprime,M);
  }
}

void SimpleSegSiev(short *s, unsigned long n, unsigned long D, unsigned long M)
/* assumes D+1 units of memory allocated at s */
/* ensures, for 0<=j<=D, that
  s[j] = 1   if n+j is coprime to all m<=M, 
  s[j] = 0   otherwise */
{
  unsigned long j, m;
  unsigned long np;
  short bn;
  
  if(M<=1) {
    for(j=0; j<=D; j++)
      s[j] = 1;
    return;
  }

  /* we start by setting s[j]=1 for n+j odd and s[j]=0 for n+j even */
  bn = n%2;
  for(j=0; j<=D-1;) {
    s[j++] = bn; s[j++]= !bn;
  }
  if(j==D) s[j]=bn;

  if(n==0) {
    if(D>=1) s[1]=0; /* since 1 is not a prime */
    if(D>=2) s[2]=1; /* since 2 is a prime */
  }

  if(n==1) {
    s[0]=0; /* since 1 is not a prime */
    if(D>=1) s[1]=1; /* since 2 is a prime */
  }

  AllocSimpleSiev(M);

  for(m=3; m<=M; m+=2)
    if(smallprime[m]) {
      np = m*((n+(m-1))/m); /* smallest multiple >=n of m */
      if(np<=m)
	  np = 2*m;             /* don't sieve out p itself! */
      if(!(np%2))
	np+=m;
      for(; np<=n+D; np += 2*m) {
	s[np-n]=0;
      }
    }
}

void SimpleSegSievP(short *s, unsigned long n, unsigned long D)
{
  unsigned long M;
  mpz_class sqt;
  mpz_class npD(n+D);

  mpz_sqrt(sqt.get_mpz_t(), npD.get_mpz_t());
  M = sqt.get_ui();

  SimpleSegSiev(s,n,D,M);
}

void SubSegSiev(short *s, unsigned long n, long D, long M)
/* assumes D+1 units of memory allocated at s */
/* ensures, for 0<=j<=D, that
  s[j] = 1   if n+j is coprime to all m<=M, 
  s[j] = 0   otherwise */
{
  mpz_class sqt;
  mpz_class Mp1(M+1);
  unsigned long Dp,Mp,j,p;
  unsigned long np;
  short *P, bn;

  /*  fprintf(stderr,"%ld %ld %ld\n",n,D,M);*/
  if(M<=1) {
    for(j=0; j<=D; j++)
      s[j] = 1;
    return;
  }
  /* we start by setting s[j]=1 for n+j odd and s[j]=0 for n+j even,
        just as in SimpleSegSiev */
  bn = n%2;
  for(j=0; j<=D-1;) {
    s[j++] = bn; s[j++]= !bn;
  }
  if(j==D) s[j]=bn;

  if(n==0) {
    s[0]==0;
    if(D>=1) s[1]=0; /* since 1 is not a prime */
    if(D>=2) s[2]=1; /* since 2 is a prime */
  }

  if(n==1) {
    s[0]=0; /* since 1 is not a prime */
    if(D>=1) s[1]=1; /* since 2 is a prime */
  }

  
  mpz_sqrt(sqt.get_mpz_t(), Mp1.get_mpz_t());
  Dp = sqt.get_ui();                     /* Dp = (int) sqrt(M) */

  /*  fprintf(stderr,"%ld\n",Dp);*/
    
  P = (short *) calloc(Dp+1,sizeof(short));

  AllocSimpleSiev(Dp);
  for(Mp=1; Mp<=M; Mp+=Dp+1) {
    /*    fprintf(stderr,"%ld\n",Mp);*/
    SimpleSegSievP(P,Mp,Dp);
    for(p=(Mp%2 ? Mp : Mp+1); p<=Mp+Dp && p<=M; p+=2)
      if(P[p-Mp]) {
	np = p*((n+p-1)/p); 	/* smallest multiple >=n of p */
	if(np<=p)
	  np = 2*p;             /* don't sieve out p itself! */
	if(!(np%2))
	  np+=p;
	for(; np<=n+D; np += 2*p)
	  s[np-n]=0;
      }
  }
  free(P);
}

void SegSiev(short *s, unsigned long n, unsigned long D)
{
  unsigned long M;
  mpz_class sqt;
  mpz_class npD(n+D);

  mpz_sqrt(sqt.get_mpz_t(), npD.get_mpz_t());
  M = sqt.get_ui();

  SubSegSiev(s,n,D,M);
}


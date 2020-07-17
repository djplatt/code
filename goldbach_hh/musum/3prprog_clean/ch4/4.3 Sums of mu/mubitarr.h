#define logb256(m) \
  (((m)&0xFFFFFFFF00000000) ? \
  (((m)&0xFFFF000000000000) ? (((m)&0xFF00000000000000) ? 8 : 7) : \
   (((m)&0xFFFFFF0000000000) ? 6 : 5)) : \
  (((m)&0xFFFFFFFFFFFF0000) ? (((m)&0xFFFFFFFFFF000000) ? 4 : 3) : \
   (((m)&0xFFFFFFFFFFFFFF00) ? 2 : ((m) ? 1: 0))));
/* computes ceil(log_256(m+1)) */


typedef unsigned long ulong;

#define bitarrl(ar,n)  (((ar)[(n)/(8*sizeof(ulong))] &		\
			 (((ulong) 1)<<((n)%(8*sizeof(ulong))))) ? 1 : 0)

void fillprimediff(unsigned char *pdiff, short *isprime, long int N)
  /* assumes N>5*/
/* pdiff[0]:=(p_2-p_1)/2 = (5-3)/2= 1, pdiff[i] = (p_{i+2}-p_{i+1})/2, 
   pdiff[k-2] = 0,
   where k is the number of primes <=N */
{
  ulong i, j, k=0, p;
  
  for(i=3, k=0; i<N; i=j) {
    for(j=i+2; j<N; j+=2)
      if(isprime[j]) 
        break;
    pdiff[k++]  = (j-i)/2;
  }
  
  if(k)
    pdiff[k-1] = 0;
}

/*The following function fills out arrays containing mu, mu^2 and other
  useful information for the first M integers, M as below.
  These small arrays will be copied time and again over larger arrays later,
  before sieving procedures start,
  so as to save us the trouble of sieving by small primes */
void ninitmudiv(ulong *imus, ulong *isqf, unsigned char *indgcd,
	        short *isprime,
		ulong minsq, ulong minpr, ulong M, ulong Mm)
/* Preconditions:

  M = \prod_{p<minsq} p  *   \prod_{p<minpr} p 
  Mm =  \prod_{p<minpr} p 
  minpr>=minsq>2
  #{primes<minpr} <=8 (i.e. minpr<=23)
  isprime[i] contains whether i is a prime for i<minpr

  imus and isqf have >= (M-1)/(8*sizeof(ulong))+1  entries
  indgcd has >= Mm entries
*/

/* Effects:

   imus, isqf will be bit-arrays; initlogpr is a hex-array
    imus, isqf have  M useful entries
   for 0<=k<Mm, 
     indgcd[k] will contain the number (b_{l-1} ... b_0) in base 2,
     where b_i=1 if p_i|k and b_i=0 otherwise
   (Here p_0=2, p_1=3,... are the primes.)
   if k is free of square divisors p^2, p<minsq:
    isqf_k = 1
    imus_k = \prod_{p|k, p<minpr} (-1)
   else
    isqf_k = 0
    imus_k has undefined behavior
*/

/* the end stubs (one ulong) of the bit- and hex-arrays are padded with 0s*/
{
  ulong j, p, p2;
  ulong bind;
  unsigned int jbit;
  unsigned char u;
  
  memset(imus,0xAA,(M-1)/8+1);
    /* 0 if = 0 or 2 mod 4, 1 if = 1 or 3 mod 4 */
  memset(isqf,0xEE,(M-1)/8+1);
    /* 0 if = 0 mod 4, 1 otherwise */
 
  for(j=0; j<M; j+=2) {
    indgcd[j] = 1;
    indgcd[j+1] = 0;
  }
     /* in other words, 1 if even, 0 if odd */
  
  for(p=3, u=0x01; p<minpr; p+=2)
    if(isprime[p]) {
      u<<=1;
      for(j=0; j<M; j+=p) {
	bind = j/(sizeof(ulong)*8);
	jbit = (j%(sizeof(ulong)*8));

	/* presumably the compiler is clever enough to do divisions
	   and moduli by powers of 2 quickly */;

	imus[bind] ^= ((ulong) 1)<<jbit;
	/* flips mu[j] */	

	if(j<Mm)
	  indgcd[j] |= u;
	/* adds 2^i to indgcd[j], if j is divisible by the ith prime */
      }
      
      if(p<minsq) {
	p2 = p*p;
	for(j=0; j<M; j+=p2) {
	  bind = j/(sizeof(ulong)*8);
	  jbit = (j%(sizeof(ulong)*8));
	  
	  isqf[bind] &= ~(((ulong) 1)<<jbit);
	  /* sets sqf[j] to 0 */
	}
      }
    }

  /*... and here we fill the end stubs with 0s */
  /*  imus[(M-1)/(8*sizeof(ulong))] &=
    (~((ulong) 0)) >> 
    ((8*sizeof(ulong))-M%(8*sizeof(ulong)));
  isqf[(M-1)/(8*sizeof(ulong))] &=
    (~((ulong) 0)) >> 
    ((8*sizeof(ulong))-M%(8*sizeof(ulong)));*/
}

void copyblock(ulong *dest, ulong *src, ulong m, ulong M)
/* Preconditions: M divides m
   as ulong arrays:
   dest has (m-1)/(8*sizeof(ulong))+1 free entries
   src  has (M-1)/(8*sizeof(ulong))+1 entries */

/* copies the M entries of the bit arrays src, isqf to dest, sqf
   m/M times */
{
  unsigned long i, j, ib0;
  unsigned int offb, nbits = 8*sizeof(ulong);
  ulong tail;
  /* the last M-((M-1)/nbits)*nbits bits */

  if(M%nbits)
    tail = src[(M-1)/nbits] & ((~ ((ulong) 0))>>(nbits - M%nbits));
  else
    tail = src[(M-1)/nbits];
  
  memset(dest,0,(m-1)/8+1);

  memcpy(dest,src,(M-1)/8);
  dest[(M-1)/nbits] = tail;
  /*    printf("%lx %lx %lu\n",tail,((~ ((ulong) 0))>>(nbits - M%nbits)),M%nbits);*/
   
  
  for(i=M; i<m; i+=M) {
    offb = i%nbits;
    ib0 = i/nbits;
    
    if(offb) {
      /*     if(i==4*M)
	     printf("%lu %u %lx \t %lx %lx\t",i,offb,src[0],dest[ib0],dest[ib0+1]);*/
      for(j=0; j<(M-1)/nbits; j++) {
	dest[ib0+j]   |= src[j]<<offb;	
	dest[ib0+j+1] |= src[j]>>(nbits-offb);
      }
      /*          if(i==4*M)
		  printf("%lx %lx\n",dest[ib0],dest[ib0+1]);*/
      dest[ib0+j] |= tail<<offb;
      if(ib0+j+1<=(m-1)/nbits) 
	dest[ib0+j+1] |= tail>>(nbits-offb);
    } else  {
      memcpy(dest+ib0,src,(M-1)/8);
      dest[ib0+(M-1)/nbits] = tail;
    }
  }
  /*              printf("%lx %lx\n",dest[M/nbits],dest[M/nbits+1]);*/
}

ulong *makegcdind(unsigned char *indgcd, short *isprime, ulong minpr)
/* returns an array ar of length 2^l, where l is # of primes <= minpr
   For n = (b_{l-1} ... b_0) in base 2,
   ar[n] will be \prod_{i=0}^{l-1} p_i^{b_i}, where p_i is the ith prime */
{
  unsigned int N, p, pind, i, u;
  ulong *gcdind;
  
  for(N=2, p=3, pind=0; p<minpr; p+=2)
    if(isprime[p])
      N <<= 1;

  gcdind = (ulong *) calloc(N, sizeof(ulong));
  
  for(i=0; i<N; i+=2) {
    gcdind[i] = 1; gcdind[i+1] = 2;
  }

  for(u=0x02, p=3, pind=0; p<minpr; p+=2)
    if(isprime[p]) {
      for(i=0; i<N; i++)
	if(i&u)
	  gcdind[i] *= p;
      u<<=1;
    }
  
  return gcdind;
}

void fillmubitblock(ulong *mus, ulong *sqf,
		 ulong *imus, ulong *isqf, unsigned char *indgcd,
		 ulong *gcdind, ulong *lor,
		 unsigned char *pdiff,
		 ulong n, ulong m,
		 ulong minsq, ulong minpr, ulong M, ulong Mm)
/* sets mus[0].. mus[m/sizeof(long)],
        sqf[0].. sqf[m/sizeof(long)] in such a way that
	the bth bit of mus[j] stores (1+mu(n+sizeof(long)*j+b))/2
            if (n+sizeof(long)*j+b) is square-free, and
        the bth bit of sqf[j] stores whether
               n+sizeof(long)*j+b is square-free */

/* assumes isprime is filled and valid up to and including sqrt(n+m-1) */
/* assumes initmu is filled and valid from 0 up to M-1,
   where M = \prod_{p<minsq} p  *   \prod_{p<minpr} p */
/* assumes Mm = \prod_{p<minpr} p */
/* assumes gcdind and indgcd are as returned by makegcdind and ninitmudiv */
/* assumes n and m are divisible by M */
/* assumes 2 < minsq <= minpr */
/* convention: mu(0) = 0 */
{
  ulong i,j,mi,di;
  ulong p, p2, maxp;
  ulong *logpr;
  ulong bind, hind, pind;
  unsigned int jbit, jhex;
  ulong band, lo;
  ulong u;
  ulong pref;
  unsigned char ig;
  ulong pd;
  
  logpr = (ulong *) calloc(m/(2*sizeof(ulong))+1,sizeof(ulong));
  
  copyblock(mus, imus, m, M);
  copyblock(sqf, isqf, m, M);

  /*    for(i=0; i<m; i+=M) {
    printf("%lu:\t",i);
    for(j=i; j<i+8; j++) {
      if(bitarrl(sqf,j))
	if(bitarrl(mus,j))
	  printf(" 1 ");
	else
	  printf("-1 ");
      else
	printf(" 0 ");
    }
    printf("\n");
    }*/

  /*  for(i=0; i<m; i+=M) {
    printf("s %lu:\t",i);
    for(j=i; j<i+8; j++) {
      if(bitarrl(sqf,j))
	printf(" 1 ");
      else
	printf(" 0 ");
    }
    printf("\n");
    printf("m %lu:\t",i);
    for(j=i; j<i+8; j++) {
      if(bitarrl(mus,j))
	printf(" 1 ");
      else
	printf(" 0 ");
    }
    printf("\n");
    }*/
      
  maxp = (ulong) sqrt(n+m-1);

  
  for(p=3,pind=0; p<minsq; p+=2*pdiff[pind++])
    ;
  /* sets p to the smallest prime p<=minsq */
  
  for(band = 0xFF, lo=1; p<=maxp; p+=2*pd) {
    for( ; p&(~band) ; lo++)
      band = (band<<8)|((ulong) 0xFF);
    /* now lo contains ceil(log_256 p);
       a prime cannot be a power of 256 after all */
    
    if(p>=minpr) {
      for(j=(n ? ((((n+p)-1)/p)*p-n) : p); j<m; j+=p) {
	/* (((n+p)-1)/p)*p is the smallest multiple of p that is >=n */
	bind = j/(sizeof(ulong)*8);
	hind = j/(sizeof(ulong)*2);
	jbit = (j%(sizeof(ulong)*8));
	jhex = (j%(sizeof(ulong)*2))*4;
	
	mus[bind] ^= (((ulong) 1)<<jbit);
	/* flips mu */

	logpr[hind] += lo<<jhex;
	  /* adds ceil(log_256 p) to logpr_j*/
	}
      }
      
    p2 = p*p;
    for(j=(n ? ((n+p2-1)/p2)*p2-n : p2); j<m; j+=p2) {
      bind = j/(sizeof(ulong)*8);
      jbit = (j%(sizeof(ulong)*8));
	    
      sqf[bind]  &= ~(((ulong) 1)<<jbit);
      /* sets j+n as non-squarefree */
    }

    if(!(pd=pdiff[pind++]))    /* last prime in our set */
      break;
  }

  /*Now comes the last step of the sieving process */
  for(i=max(n,1), mi=max(n,1)%Mm, bind=hind=0,
        jhex=(n ? 0 : 4), u = (n ? 1 : 2); i<n+m; i++) {
    if(sqf[bind] & u) {     /* if i is square-free... */
      if(ig=indgcd[mi]) {  /* if i has prime factors < minpr... */
	  if((pref=gcdind[ig])<16) 
	 /* pref=gcdind[ig] is the product of the prime factors <minpr of i */
	    {
	    lo = ((i>(pref<<(4*8))) ?
		  ((i>(pref<<(6*8))) ? ((i>(pref<<(7*8))) ? 8 : 7) :
		   ((i>(pref<<(5*8))) ? 6 : 5)) :
		  ((i>(pref<<(2*8))) ? ((i>(pref<<(3*8))) ? 4 : 3) :
		   ((i>(pref<<(1*8))) ? 2 : (i>pref ? 1 : 0))));
	    /* lo = ceil(log_256(i/pref)) */
	  } else {
	    di = i-1;
	    lo = logb256(di);  
	    lo -= lor[ig];   /* lo = ceil(log_256(i)) - ceil(log_256(pref)) */
	  }
	} else {
	  di = i-1;
	  lo = logb256(di);  /* lo = ceil(log_256(i)) */
	}

      /*      if(i==180183)
	if(((logpr[hind]>>jhex)&((ulong) 0x0F))<lo)
	  if(mus[bind]  ^ u)     
	    printf("FLIP to 1!\n");
	  else
	    printf("FLIP to -1!\n");
	else
	printf("NO FLIP!\n");*/
      
      if(((logpr[hind]>>jhex)&((ulong) 0x0F))<lo)
	/* iff the sum of ceil(log_256(p)) for all primes>=minpr in our set
            is < lo, then there is a missing prime factor and so... */
	mus[bind]  ^= u;     /* flip mu[i] */

      /* in that iff, "if" is clear:
          we have \sum_{p>=minpr in our set} ceil(log_256(p)) 
          <ceil(log_256(i/pref)) or < ceil(log_256(i)) - ceil(log_256(pref)),
	  and so i has to be larger than the product of pref and those primes.
        as for "only if": 
	  We are using the assumption that minpr>=17, and also that
	  we have all the primes up to sqrt(n+M-1), inclusive.
	  cei((log_256(m))<=2 log_256(m) for any m>=16.
	  Hence \sum_{p>=minpr in our set} ceil(log_256(p)) <=
	       2 \sum_{p>=minpr in our set} log_256(p) <=
	       (log_256((product of such primes)^2)), and
	       that is less than log_256(i) if i has a prime factor 
	       outside the set (and thus > sqrt(i))
          The same goes for 
	  \sum_{p>=minpr in our set} ceil(log_256(p)) + ceil(log_256(pref))
            if pref>=16 */
    }

    /* now we just increment all the indices */
    u<<=1; jhex+=4; mi++;
    
    if(mi==Mm) mi = 0;
      
    if(jhex==(8*sizeof(ulong))) {
      jhex=0; hind++;
    }
    
    if(!u) {
      bind++; u=1;
    }  
  }
  
  free(logpr);
}

void initmubitdat(ulong N, ulong minsq, ulong minpr,
		ulong **imus, ulong **isqf, unsigned char **indgcd,
		ulong **gcdind, ulong **lor,
		unsigned char **pdiff, 
		ulong *M, ulong *Mm)
/* requires: minsq>=3, minpr>=17 both prime*/
/* allocates and initializes data necessary for working efficiently
   with bitarrays containing mu */
{
  short *isprime;
  ulong p, K, i,j,m;
  ulong nprimes, maxp;

  isprime = (short *) calloc(minpr+1,sizeof(short));
  fillisprime(isprime,minpr+1); 

  for(*M = 1, p=2; p<minsq; p++)
    if(isprime[p])
      *M *= p;
  
  for(*Mm = 1, p=2, K=1; p<minpr; p++)
    if(isprime[p]) {
      *M *= p; *Mm *= p; K<<=1; 
    }

    
  *imus = (ulong *) calloc(*M/(8*sizeof(ulong))+1,sizeof(ulong));
  *isqf = (ulong *) calloc(*M/(8*sizeof(ulong))+1,sizeof(ulong));
  *indgcd = (unsigned char *) calloc(*Mm,sizeof(ulong));

  ninitmudiv(*imus, *isqf, *indgcd, isprime, minsq, minpr, *M, *Mm);
  *gcdind = makegcdind(*indgcd,isprime,minpr);
  
  free(isprime);

  maxp = (ulong) (-sqrt((int_double) N+1).right);

  isprime = (short *) calloc(maxp+1,sizeof(short));
  fillisprime(isprime,maxp+1); 
  if(maxp>=2) nprimes=1; else nprimes=0;
  for(i=3, nprimes=0; i<=maxp; i+=2)
    if(isprime[i])
      nprimes++;
  *pdiff = (unsigned char *) calloc(nprimes,sizeof(unsigned char));
  fillprimediff(*pdiff,isprime,maxp+1);
  free(isprime);

 *lor = (ulong *) calloc(K,sizeof(ulong));

  for(i=0; i<K; i++) {
    j = (*gcdind)[i]-1;
    (*lor)[i] = logb256(j);
    /* so lor[i] stores ceil(log_256(gcdind[i])) */
  }
}

void freemubitdat(ulong *imus, ulong *isqf, unsigned char *indgcd,
		  ulong *gcdind, ulong *lor,
		  unsigned char *pdiff)
{
  free(imus); free(isqf); free(indgcd); free(gcdind); free(lor); free(pdiff);
}

void fillmufrombitarrs(short *arr, ulong *mus, ulong *sqf, int L)
{
  int i;

  for(i=0; i<L; i++)
    if(bitarrl(sqf,i))
      if(bitarrl(mus,i))
	*arr++ = 1;
      else
	*arr++ = -1;
    else
      *arr++ = 0;
}


void printBits(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char bit;
    int i, j;

    for (i=0;i<size;i++)
    {
        for (j=0;j<8;j++)
        {
            bit = (b[i] >> j) & 1;
            printf("%u", bit);
        }
    }
}

void printmubitdat(ulong N, ulong minsq, ulong minpr,
		   ulong *imus, ulong *isqf,
		unsigned char *pdiff, 
		ulong M)
{
  ulong p,m;
  
  printf("2\t");
  for(p=3; *pdiff; p+=2*(*pdiff++))
    printf("%lu\t",p);
  printf("%lu\n",p); 
  for(m=0; m<(M-1)/(8*sizeof(ulong))+1; m++) {
    printBits(sizeof(ulong),&(imus[m])); printf("\t");
    printBits(sizeof(ulong),&(isqf[m])); printf("\n");
    }
}

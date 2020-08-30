
void fillmdot(int_double *mdot, short *isprime,
	      long int d, long int y, short int *mu, long int *sigma)
/* stores mdot_{d}(n) in h[n] for all 1<=n<=y */
/* assumes mu[n] and sigma[n] are filled correctly for 1<=n<=y*/
{
  long i;
  short int *cop;

  cop = (short int *) calloc(y+1,sizeof(short int));
  fillcopblock(cop,isprime,0,y+1,d);
  mdot[0]=0.0;
  for(i=1; i<=y; i++) {
    if(cop[i] && mu[i])
      mdot[i] = mdot[i-1] + ((int_double) mu[i])/sigma[i];
    else
      mdot[i] = mdot[i-1];
  }
  free(cop);
}

void fillmtilde(int_double *mtilde, int_double *mdot, long int d, long int y, short int *mu, long int *sigma, long int sigmad)
/* stores mtilde_{d}(n)-zeta(2) sigma(d)/d in mtilde[n] for all 1<=n<=y */
/* assumes mu[n], sigma[n] and mdot[n] (= mdot_{d}(n)) are filled correctly for 1<=n<=y*/
{
  long i;

  mtilde[0]=mtilde[1]= - (((int_double) sigmad)/d)*d_pi*d_pi/6.0;
  for(i=2; i<=y; i++)
    mtilde[i] = mtilde[i-1] + mdot[i-1]*(-log1p(-1/((int_double) i)));
}

void fillmdmtarray(int_double ***mdv, int_double ***mtv,
		   int v, int M,
		   short *isprime, short int *mu, long int *sigma)
/* Allocates arrays *mdv, *mtv of length M+1
   and, for each 1<=d<=M, arrays (*mdv)[d], (*mtv)[d] of length M/d+1 */
/* Stores mtilde_{d v}(n) - zeta(2) sigma(d v)/d v in (*mdv)[d][n] 
      and mdot_{d v}(n) in (*mdot)[d][n]  
   for 1<=d<=M  and 1<=n<=M/d */
/* assumes isprime[n], mu[n], sigma[n] correctly initialized for 1<=n<=M */
{
  int d;
  
  (*mdv)=(int_double **) calloc(M+1,sizeof(int_double *));
  (*mtv)=(int_double **) calloc(M+1,sizeof(int_double *));
  for(d=1; d<=M; d++) {
    (*mdv)[d] = (int_double *) calloc(M/d+1,sizeof(int_double));
    (*mtv)[d] = (int_double *) calloc(M/d+1,sizeof(int_double));
    fillmdot((*mdv)[d],isprime,v*d,M/d,mu,sigma);
    fillmtilde((*mtv)[d],(*mdv)[d],v*d,M/d,mu,sigma,sigma[d]*sigma[v]);
  }
}
		 
void fillm(int_double *m, int v, long int y, short int *mu)
/* stores m_v(n) in m[n] for all 1<=n<=y (v=1 or v=2)*/
/* assumes mu[n] is filled correctly for 1<=n<=y*/
{
  long i;
  
  m[0]=0.0;
  for(i=1; i<=y; i+=v) {
    if(mu[i])
      m[i] = m[i-1] + ((int_double) mu[i])/i;
    else
      m[i] = m[i-1];
    if(v==2 && i<y)
      m[i+1] = m[i];
  }
}

void fillmch(int_double *mch, int_double *m, int v, long int y,
             short int *mu)
/* stores mcheck_{v}(n)- v/phi(v) in mch[n] for all 1<=n<=y (v=1 or v=2)*/
/* assumes mu[n], sigma[n] and m[n] (= m_{v}(n)) 
   are filled correctly for 1<=n<=y. Oh, and y>=1*/
{
  long i;

  mch[0]= mch[1]= -v;
  for(i=2; i<=y; i+=v) 
    mch[i] = mch[i-1] + m[i-1]*(-log1p(-1/((int_double) i)));
}


void siep(primfact *factar, long d, primfact **newf)
/* select the primes within a list that do *not* divide d */
{
  primfact *fa;

  for(*newf = NULL; factar; factar = factar->next)
    if(d%(factar->p)) {
      fa = (primfact *) malloc(sizeof(primfact));
      fa->p = factar->p; fa->r = factar->r;
      fa->next = *newf;
      *newf = fa;
    }
}

void purg(primfact *factar)
{
  if(factar) {
    purg(factar->next);
    free(factar);
  }
}

void fillk(int_double *k0, int_double *k1, int_double *k2,
	   int_double *k0init, int_double *k1init, int_double *k2init,
	       int v, long n0, long n1, long Sqt,
	       int_double *md2, int_double *mt2,
	       short int *isprime, short int *mu,
	   long int *sigma, int_double **mdv, int_double **mtv)
/* fills 
k0[0], k0[1],... k0[n1-(n0+1)] 
k1[0], k1[1],... k1[n1-(n0+1)] 
k2[0], k2[1],... k2[n1-(n0+1)] 
with
   k_{0,v}(n0+1),...,k_{0,v}(n1)   
   k_{1,v}(n0+1),...,k_{1,v}(n1) 
   k_{2,v}(n0+1),...,k_{2,v}(n1) */
/* Preconditions:
    v = 1 or 2
    Sqt = floor of sqrt(n1)
    isprime[1],...,isprime[Sqt] filled correctly 
    mu[1],...,mu[Sqt] filled correctly ( = mu(1), mu(2),...)
    sigma[1],...,sigma[Sqt] filled correctly (= sigma(1), sigma(2),...)
    if n0>0:
      *k0init, *k1init, *k2init set to k_{0,v}(n0), k_{1,v}(n0), k_{2,v}(n0)
      md2(d) set to \dot{m}(floor(n0/d))
      mt2(d) set to \tilde{m}(floor(n0/d))
*/
  
{
  primfact **factar;
  short int *mu2;
  short int *cop;
  long n, i;
  int_double k1sum,k2sum,k1ss,k2ss,L;
  long int sigma2;
  long int *divo, *d, *sigo, *s, *ldiv, *lsig;
  primfact *fa, *facy;
  long int *l, *ls, r;
  int_double constter, mfact;

  factar = (primfact **) calloc(n1-n0+1,sizeof(primfact *));
  mu2 = (short int *) calloc(n1-n0+1,sizeof(short int));
  cop = (short int *) calloc(n1-n0+1,sizeof(short int));
     
  fillfactblock(factar, isprime, n0+1, n1-n0+1);
  fillmublock(mu2, isprime, n0+1, n1-n0+1);
  fillcopblock(cop, isprime, n0+1, n1-n0+1, v);

  if(v==1) {
    constter = d_pi*d_pi/6.0;
    mfact = 2*d_pi*d_pi/6.0;
  }  else {
    constter = (3.0*d_pi*d_pi)/6.0;
    mfact = (3.0*d_pi*d_pi)/6.0;
  } /*constter = zeta(2) sigma(v)/phi(v), mfact = 2 zeta(2) sigma(v)/v */
  

   
  if(!n0) {
    *k2init = k2[0] = 1.0;
    *k0init = k0[0] =  constter;   
    *k1init = k1[0] = -mfact;
    md2[1] = 1.0;
    for(i=v+1; i<=Sqt; i+=v)     
      md2[i] = 0;
    for(i=1; i<=Sqt; i+=v)     
      mt2[i] = -(d_pi*d_pi/6.0)*((int_double) (sigma[i]*sigma[v]))/(i*v);
  }

  for(n=max(2,n0+1); n<=n1; n++)     {
    divsiglist(factar[n-(n0+1)],&divo,&sigo);

    for(d=divo; *d; d++) 
      if(*d<=Sqt) 
	if(*d < n)
    	  mt2[*d] += md2[*d]*log1p(((int_double) *d)/(n-(*d)));
    /* \tilde{m}_{dv}(n/d) = \tilde{m}_{dv}(n/d-1) 
                            + \tilde{m}_{dv}*log(n/(n-d)) */

    
    if(mu2[n-(n0+1)] && cop[n-(n0+1)]) {
      for(fa=factar[n-(n0+1)], sigma2 = n; fa; fa=fa->next)
	sigma2 += sigma2/(fa->p);
      /* sigma2 = n * \prod_{p|n} (1+1/p) = sigma(n) (since n is squarefree)*/
 
      k1sum = k2sum = 0.0;
      for(d=divo, s=sigo; *d; d++, s++) {	
	if(*d<=Sqt) {
	  k2sum += md2[*d]/(*s);
	  /* add \dot{m}_{d v}(n/d - 1)/sigma(d) */
	  k1sum += mt2[*d]/(*s);
	  /* add \tilde{m}_{d v}(n/d)/sigma(d) */
	}	     
	if((*d)*Sqt<n) {
	  k1ss=k2ss=0.0;
	  siep(factar[n-(n0+1)],*d,&facy);
	  divsiglist(facy, &ldiv, &lsig);
	  /*facy lists the prime divisors of n/d */
	  /*ldiv lists the divisors l of n/d, and lsig lists sigma(l) */
	  
	  for(l=ldiv, ls=lsig; *l; l++, ls++) {
	    if(*l<*d) {
	      r = (*d)/(*l);
	      k1ss+=(mtv[*l][r] +
		     mdv[*l][r] *
		     log(((int_double) *d)/(r*(*l))))/(*ls);
	      k2ss+=mdv[*l][(*d-1)/(*l)]/(*ls);

	      /* k1ss += (1/sigma(l))*(\tilde{m}_{l v}(d/l)
                                        - (sigma(lv)/(lv)) zeta(2))         
		 k2ss += (1/sigma(l))*\dot{m}_{l v}((d-1)/l) */
	    } else {
	      k1ss-=mfact*0.5/(*l);
	      /* k1ss -= zeta(2) sigma(v)/l v */
	    }
	  }
	  
	  k1sum+=k1ss/((int_double) (sigma2/(*s)));
	  k2sum+=k2ss/((int_double) (sigma2/(*s)));
	  /* here (sigma2/(*s)) = sigma(n)/sigma(d) */ 
	  purg(facy); free(ldiv); free(lsig);
	}
      }


      /*now we should update \dot{m}_{dv}(n/d-1) to \dot{m}_{dv}(n/d) */ 
      for(d=divo, s=sigo; *d; d++, s++) 
	if(*d<=Sqt)
	    /* if n is square-free and (n,v)=1, then (n/d,d)=1 and (n/d,v)=1;
	       if n is not square-free, then either n/d is not
	       square-free, or (n/d,d)\ne 1, or d is not square-free,
	       and so md2[*d] need not be updated.
	       if n is square-free and (n,v)\ne 1, then either
                 (n/d,v)\ne 1 or (d,v)\ne 1,
		 and so md2[*d] need not be updated.*/
	    md2[*d] += (mu2[n-(n0+1)]*mu[*d])/((int_double) (sigma2/(*s)));
      /* adds mu(n/d)/sigma(n/d) */
	
	L = -log1p(-((int_double) 1.0)/n);
	*k0init =
	  k0[n-(n0+1)] = *k0init + (*k1init)*L + (*k2init)*L*L;
	*k1init =
	  k1[n-(n0+1)] = *k1init + 2*(*k2init)*L
	  + 2*k1sum*mu2[n-(n0+1)]/((int_double) sigma2);
	*k2init =
	  k2[n-(n0+1)] = *k2init + k2sum*2*mu2[n-(n0+1)]/((int_double) sigma2);
    } else {
	L = -log1p(-((int_double) 1.0)/n);
	*k0init =
	  k0[n-(n0+1)] = *k0init + (*k1init)*L + (*k2init)*L*L;
	*k1init = k1[n-(n0+1)] = *k1init+2*(*k2init)*L;
	k2[n-(n0+1)] = *k2init;
    }
       free(divo); free(sigo);

  }
  
  for(i=0; i<=n1-n0; i++)
    purg(factar[i]);
  
  free(factar); free(mu2); free(cop);
} 

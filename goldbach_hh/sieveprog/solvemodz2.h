
typedef struct intlist intlist;

struct intlist {
  long a, q, end;   /* a, a+q, a+2q... while terms <=end */
  intlist *next;
};

inline long longceil(mpq_class x)
{
  mpz_class res;
  
  mpz_cdiv_q(res.get_mpz_t(), x.get_num_mpz_t(), x.get_den_mpz_t());

  return res.get_si();
}

inline long longfloor(mpq_class x)
{
  mpz_class res;
    
  mpz_fdiv_q(res.get_mpz_t(), x.get_num_mpz_t(), x.get_den_mpz_t());

  return res.get_si();
}

intlist *inclist(long a, long q, long end, intlist *l)
{
  intlist *head;


  head = (intlist *) malloc(sizeof(intlist));
  if(!head)
    exit(1);
  head->a = a;   head->q = q;   head->end = end;
  head->next = l; 


  return head;
}

intlist *intintertolist(long Imin, long Iplu)
  /* returns the list of integers in the closed interval [Imin, Iplu] */
{
  intlist *head;

  if(Imin>Iplu)
    return NULL;
  else 
    return inclist(Imin,1,Iplu,NULL);
}

void killlist(intlist *l)
{
  if(l) {
    killlist(l->next);
    free(l);
  }
}

intlist *lastlist(intlist *l)
{
  intlist *last;

  if(!l)
    return NULL;

  for(last=l; last->next; last = last->next)
    ;
  return last;
}

intlist *appendlist(intlist *l1, intlist *l2)
{
  intlist *l;
  
  if(!l1)
    return l2;
  else {
    l=lastlist(l1);
    l->next = l2;
    return l1;
  }
}

/*intlist *mergelist(intlist *l1, intlist *l2)*/
  /* merge two lists already in increasing order */
/*{
  if(!l1)
    return l2;
  if(!l2)
    return l1;
  
  if(l1->it<=l2->it) {
    if(l1->it==l2->it) {
      l1->next = mergelist(l1->next,l2->next);
      free(l2);
    } else
      l1->next = mergelist(l1->next,l2);
    return l1;
  } else {
    l2->next = mergelist(l2->next,l1);
    return l2;
  } 
}
*/

intlist *intertolist(mpq_class Imin, mpq_class Iplu)
/* returns the list of integers in the closed interval [Imin, Iplu] */
{
  /*  gmp_fprintf(stderr,"from %Qd to %Qd\n",Imin.get_mpq_t(),Iplu.get_mpq_t());
   */
   if(cmp(Imin,Iplu)<=0)
     return intintertolist(longceil(Imin),
			   longfloor(Iplu));
   else {
     /*     fprintf(stderr,"forget it\n");*/
     return NULL;
   }
}

/*void printlist(intlist *l, long i)
{
  if(l) {
     printf("%ld\t",l->it);
     printlist(l->next,i+1);
  } else
    printf("\n");
    }*/

int loud=0;

intlist *SolveInR(mpq_class alph1, mpq_class alph0,
		  mpq_class Imin, mpq_class Iplu,
		  mpq_class Jmin, mpq_class Jplu)
{
  mpq_class xmin, xplu;

  if(loud)
    gmp_fprintf(stderr,"SolveInR: %Qd %Qd %Qd %Qd %Qd %Qd\n",
		alph1.get_mpq_t(),alph0.get_mpq_t(),
	      Imin.get_mpq_t(),Iplu.get_mpq_t(),
	      Jmin.get_mpq_t(),Jplu.get_mpq_t());

  if(sgn(alph1)==0) {
    if(cmp(Jmin,alph0)<=0 && cmp(alph0,Jplu)<=0) {
	return intertolist(Imin,Iplu);
      }
    else {
      return NULL;
    }
  } else {
    if(sgn(alph1)>=0) {
      xmin = Jmin - alph0;  xplu = Jplu - alph0;
    } else {
       xplu = Jmin - alph0;  xmin = Jplu - alph0;
    }
    xmin /= alph1;  xplu /= alph1;
    
    return intertolist((cmp(Imin, xmin)>0 ? Imin : xmin),
		       (cmp(Iplu, xplu)<0 ? Iplu: xplu));
  }
}

intlist *SmallSolveModZ(mpq_class alph1, mpq_class alph0,
			mpq_class Imin, mpq_class Iplu,
			mpq_class eta)
{
  mpq_class Lmin, Lplu, xmin, xplu;
  mpz_class nmin, nplu;
  intlist *r, *l;

  if(loud)
    gmp_fprintf(stderr,"SmallSolveModZ: %Qd %Qd %Qd %Qd %Qd\n",
		alph1.get_mpq_t(),alph0.get_mpq_t(),
		Imin.get_mpq_t(),Iplu.get_mpq_t(),eta.get_mpq_t()); 

  if(sgn(alph1)>=0) {
    Lmin = alph1*Imin+alph0;    Lplu = alph1*Iplu+alph0;
  } else {
    Lplu = alph1*Imin+alph0;    Lmin = alph1*Iplu+alph0;
  }

  if(loud)
    gmp_fprintf(stderr,"%.10g %.10g\n",
		Lmin.get_d(),Lplu.get_d());
  Lmin += eta; Lplu +=eta;
  mpz_fdiv_q(nmin.get_mpz_t(),
	     Lmin.get_num_mpz_t(),
	     Lmin.get_den_mpz_t());
  mpz_fdiv_q(nplu.get_mpz_t(),
	     Lplu.get_num_mpz_t(),
	     Lplu.get_den_mpz_t());
  if(loud)
    gmp_fprintf(stderr,"%Zd %Zd\n",
		nmin.get_mpz_t(),nplu.get_mpz_t());
  r = SolveInR(alph1,alph0,Imin,Iplu,nmin-eta,nmin+eta);
  if(cmp(nplu, nmin)>0) {
    l = lastlist(r);
    if(l) 
      l->next = SolveInR(alph1,alph0,Imin,Iplu,nplu-eta,nplu+eta);
    else
            r = SolveInR(alph1,alph0,Imin,Iplu,nplu-eta,nplu+eta);
  }
  
  return r;
}



intlist *SolveModZ(rat alpha1, rat alpha0,
		   long R, mpq_class eta)
{
  long a, ainv, q, j, rho, prerho, k, c, rl;
  intlist *r, *r2, *s, *head;
  mpq_class alph1(alpha1.num, alpha1.den);
  mpq_class alph0(alpha0.num, alpha0.den);
  mpq_class half(1,2);
  mpq_class delta;
  mpz_class res, resrat;
  
  diophappr(alpha1, 2*R, &a, &ainv, &q);
  if(loud)
    fprintf(stderr,"a = %ld, ainv = %ld, q = %ld\n",a,ainv,q); 
  mpz_fdiv_r(res.get_mpz_t(),
	     alph0.get_num_mpz_t(), alph0.get_den_mpz_t());
  mpq_class alph0frac(res, alph0.get_den());

  
  resrat = alph0frac*q+half;
  c = longfloor(resrat);
  delta = q*alph1 - a;
  k = longfloor(eta*q);

  if(loud) 
    fprintf(stderr,"c = %ld, k = %ld, delta = %g\n",c,k,
	    delta.get_d());
  r = NULL; r2 = NULL;
    
  for(j=-k-1; j<=k+1; j++) {
    prerho = -ainv*(c+j); rho = mod(prerho,q);
    if(loud)
      fprintf(stderr,"j = %ld\t prerrho = %ld\t rho = %ld\n",
	      j, prerho, rho);
    /*        fprintf(stderr,"j = %ld, k = %ld, c = %ld, rho = %ld\n",j,k,c,rho);*/
    if(j>=-k+1 && j<=k-1) {
      rl = div(-R-rho,q)*q+rho; if(rl<-R) rl+=q;
      r = inclist(rl,q,R,r);
    } else {
      mpq_class Imin(-R-rho, q);
      mpq_class Iplu(R-rho,q);
      Imin.canonicalize(); Iplu.canonicalize();
      
      s = SmallSolveModZ(delta, rho*alph1+alph0, Imin, Iplu, eta);
      for(head=s; head; head = head->next)
	r2=inclist(q*(head->a)+rho,q*(head->q),q*(head->end)+rho,r2);
      if(s)
	killlist(s);
    }
  }
  r =   appendlist(r,r2);
  return r;
}

/*void printxylist(intlist *list, mpq_class alph1, mpq_class alph0)
{
  long r;
  mpq_class res;
  mpz_class red;
  double resd;
  
  for(; list; list = list->next) {
    r = list->it;
    res =alph1*r+alph0;
    
    mpz_fdiv_r(red.get_mpz_t(), res.get_num_mpz_t(), res.get_den_mpz_t());
    res.get_num() = red;
    if((res.get_den()).get_si())
      resd = ((double) (res.get_num()).get_si())/
	(res.get_den()).get_si();
    else
	  fprintf(stderr,"waaaaah\n");
    gmp_printf("%ld %g\n",r,resd);
  }
}
*/

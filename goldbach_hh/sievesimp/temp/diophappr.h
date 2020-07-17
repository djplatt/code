typedef struct {
  long num, den;
} rat;

#define mod(x,y) ((x)>=0 ? ((x)%(y)) : ((y)-1) - ((-(x)-1)%(y)))
#define div(x,y) ((x)>=0 ? ((x)/(y)) : -1-(-(x)-1)/(y))

void diophappr(rat alph, long Q, long *a, long *ainv, long *qout)
/* precondition: alph.num>=0, alph.den>=0 */
/* constructs approximation a/q, q<=Q, to alph */
/* sets ainv to a^{-1} mo q */
{
  long b, p, q, pmin, qmin, pplus, qplus, nummodb, flip;
  int s;

  b = alph.num/alph.den;
  p = b; q = 1; pmin = 1; qmin = 0; s = 1;
  while(q<=Q) {
    /*    printf("%ld %ld\n",p,q);*/
    nummodb = alph.num%alph.den;
    if(nummodb==0) {
      flip = (s==1 ? -qmin : qmin);
      *a = p; *ainv = mod(flip,q); *qout = q; return;
    }
    alph.num = alph.den; alph.den = nummodb;
    b = alph.num/alph.den;

    pplus = b*p+pmin; qplus = b*q + qmin;

    pmin = p; qmin = q; p = pplus; q = qplus; s = - s;
  }

  flip = (s==1 ? q : - q);
  *a = pmin; *ainv = mod(flip,qmin); *qout = qmin;
}


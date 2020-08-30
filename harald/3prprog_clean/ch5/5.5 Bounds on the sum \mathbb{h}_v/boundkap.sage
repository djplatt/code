cei=(lambda x,d : N((((RBF(10^d)*RBF(x)).upper()).ceil())/10^d))
import re
def pf(x):
    return sage.misc.latex.LatexExpr(re.sub(r'(\d)\.?0+($|\ )',r'\1\2',latex(RR(x))).replace(r'\times',r'\cdot'))

from sage.symbolic.expression_conversions import polynomial
from sympy import sieve

def sumposcoeff(P):
  """returns the sum of positive coefficients of P"""
  Ppol = polynomial(P, base_ring=QQ)
  sum = 0
  for mon in Ppol.monomials():
    x=Ppol.monomial_coefficient(mon)
    if x>=0:
      sum+=x
  return sum

var('x','y1','y2')

def RP(k):
    P = (1-y1+x)*(1-y2+x)+y1*y2; M=y1;
    for i in [1..k]:
       M *= y2; P *= (1-M)
    P -= (1-y1+x)*(1-y2+x)
    return simplify(expand(P))

def kap(bet1,bet2,k,N):
  Rk(x,y1,y2) = RP(k)
  fk(x) = (1 + Rk(1/x,1/x^bet1,1/x^bet2)/((1-1/x^bet1+1/x)*(1-1/x^bet2+1/x)))  
  prod = 1; primes = list(sieve.primerange(1, N));
  for p in primes:
    prod *= RBF(fk(p))
  for j in [1..k]:
    prod *= RBF(zeta(bet1+j*bet2))
  alph = min(2*bet1+bet2,bet1+(k+1)*bet2); Ck=sumposcoeff(Rk)
  ktail = RBF(exp((Ck*(N-1)^(-alph+1)/(2*(alph-1)))/((1-N^(-bet1)+1/N)*(1-N^(-bet2)+1/N))))
  ktail = RBF(RIF(1,ktail.upper()))
  return prod*ktail
  
def kap2(bet1,bet2,N):
  R332(x,y1,y2) = simplify(expand((1-y1*y2)*(1-y1*y2*y2)*(1-y1*y2*y2*y2)*(1-y1*y1*y2)*(1-y1*y1*y1*y2)*((1-y1+x)*(1-y2+x)+y1*y2)-(1-y1+x)*(1-y2+x)))  
  f332(x) = (1 + R332(1/x,1/x^bet1,1/x^bet2)/((1-1/x^bet1+1/x)*(1-1/x^bet2+1/x)))  
  prod = 1; primes = list(sieve.primerange(1, N));
  for p in primes:
    prod *= RBF(f332(p))
  prod *= RBF(zeta(bet1+bet2)*zeta(2*bet1+bet2)*zeta(bet1+2*bet2)*zeta(bet1+3*bet2)*zeta(3*bet1+bet2))
  alph = min(2+bet1+bet2,1+2*bet1+bet2,1+bet1+2*bet2,4*bet1+bet2,bet1+4*bet2); C332=sumposcoeff(R332)
  ktail = RBF(exp((C332*(N-1)^(-alph+1)/(2*(alph-1)))/((1-N^(-bet1)+1/N)*(1-N^(-bet2)+1/N))))
  return (prod*ktail).upper()
  
kv = kap(1,1-1/log(10^12),1,20000)
print("kap(1,1-1/log(10^12)) is at most "+pf(cei(kv,6)))

kv = kap(1-1/log(10^12),1-1/log(10^12),1,20000)
print("kap(1-1/log(10^12),1-1/log(10^12)) is at most "+pf(cei(kv,6)))

kv = kap(1,1/2,2,40000)
print("kap(1,1/2) is at most "+pf(cei(kv,6)))

kv = kap(1-1/log(10^12),1/2,2,80000)
print("kap(1-1/log(10^12),1/2) is at most "+pf(cei(kv,6)))

kv = kap(1,1/4,4,2500000)
print("kap(1,1/4) is at most "+pf(cei(kv,6)))

kv = kap2(8/15,8/15,500000)
print("kap(8/15,8/15) is at most "+pf(cei(kv,5)))

kv = kap2(4/7,1/2,500000)
print("kap(4/7,1/2) is at most "+pf(cei(kv,5)))
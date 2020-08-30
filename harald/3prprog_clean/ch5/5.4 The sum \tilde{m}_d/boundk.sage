cei=(lambda x,d : N((((RBF(10^d)*RBF(x)).upper()).ceil())/10^d))
import re
def pf(x):
    return sage.misc.latex.LatexExpr(re.sub(r'(\d)\.?0+($|\ )',r'\1\2',latex(RR(x))).replace(r'\times',r'\cdot'))

from sage.symbolic.expression_conversions import polynomial
from sympy import sieve

var('x','y')

def sumposcoeff(P):
  """returns the sum of positive coefficients of P"""
  Ppol = polynomial(P, base_ring=QQ)
  sum = 0
  for mon in Ppol.monomials():
    x=Ppol.monomial_coefficient(mon)
    if x>=0:
      sum+=x
  return sum

def k(bet,N):
  P2(x,y) = simplify(expand((1+x-y)*(1-x*y)*(1-x*y*y) - (1+x)*(1-y)))
  f2(x) = (1 + P2(1/x,1/x^bet)/((1+1/x)*(1-1/x^bet)))  
  prod = 1; primes = list(sieve.primerange(1, N));
  for p in primes:
      prod *= RBF(f2(p))
  prod *= RBF(zeta(1+bet))*RBF(zeta(1+2*bet))
  alph = min(3*bet+1,bet+2); C2=sumposcoeff(P2)
  ktail = RBF(exp((C2*(N-1)^(-alph+1)/(2*(alph-1)))/((1+1/N)*(1-N^(-bet)))))
  ktail = RBF(RIF(1,ktail.upper()))
  return prod*ktail

bet = 1-1/log(10^12)
kv = k(bet,3000);
print("k(1-1/log(10^12)) is at most "+pf(cei(kv.upper(),6)))
    
for bet in [4/7,8/15,1/2]:
    kv = k(bet,30000);
    print("k("+repr(bet)+") is at most "+pf(cei(kv.upper(),6)))

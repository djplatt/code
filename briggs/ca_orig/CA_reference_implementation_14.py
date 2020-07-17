#!/usr/bin/env python
# Keith Briggs 2015-06-14 12:36 add sigma calculation ./CA_reference_implementation_14.py | cut -f1,6 | p -y 0 0.2
# Keith Briggs 2012-03-08 15:41 still trying a proper tail of ones!
# Keith Briggs 2012-03-07 16:43 try a proper tail of ones
# Keith Briggs 2012-03-06 15:55 Primegen
# Keith Briggs 2012-02-28 11:38 tail of ones (printed, but not stored)
# Keith Briggs 2012-02-27 17:19 cleanup
# Keith Briggs 2012-02-27 16:35 records
# Keith Briggs 2012-02-27 16:03 power sums
# Keith Briggs 2012-02-27 14:00 cleanup
# Keith Briggs 2012-02-27 13:28 keep list of epsilons
# A uses (p^3-1)/(p^3-p)=(p-1)*(p^2+p+1)/p/(p+1)/(p-1)=(1+p+1/p)/(1+p)=1+1/p/(p+1)
# B,C use (r^(a+1)-1)/(r^a-1)/r=1+1/r/(1+r+r^2+...+r^(a-1))

#from copy import copy
from math import log1p,log,floor,exp
from array import array
from sys import exit,path,stderr
#path.append('/home/kbriggs/Lagarias')
from primes import p as primes

class Primegen:
  def __init__(s):
    s.i=0
  def next(s):
    p=primes[s.i]
    s.i+=1
    return p

class P: # records for one prime
  def __init__(s,p,exponent,epsilon,power_sum):
    s.p=p
    s.exponent=exponent
    s.epsilon=epsilon
    s.power_sum=power_sum
  def update(s,exponent,epsilon,power_sum):
    s.exponent=exponent
    s.epsilon=epsilon
    s.power_sum=power_sum

class Ps: # records for all primes
  def __init__(s):
    s.ps=[P(2,1,log(7.0/6.0)/log(2),3.0),]
    s.ones=0 # first one is not counted
    s.q=P(3,0,log(4.0/3.0)/log(3),1.0) # next prime after tail of ones
  def __len__(s):
    return 1+len(s.ps) #+s.ones # count q FIXME counts ones
  def __getitem__(s,i):
    n=len(s.ps)
    if i<n: return s.ps[i]
    return s.q
  def new_q(s,q):
    s.ps.append(s.q) # move old q to list
    s.q=q
  def update(s,k,exponent,epsilon,power_sum):
    #if exponent==1: # !!!!!!!!!!!!!!!!!!!!!!!!
    #  s.ones+=1
    #  return 
    if k<len(s.ps):
      s.ps[k].update(exponent,epsilon,power_sum)
    else:
      s.q.update(exponent,epsilon,power_sum) # must be q
  def get_changeable_places(s):
    last,r=s.ps[0].exponent+1,[]
    for i,q in enumerate(s.ps):
      pe=q.exponent
      if pe<last:
        r.append(i)
        last=pe
    r.append(len(s.ps)) # allow for q
    return r
  def show_exponent_list(s):
    r='['
    ones=0
    for p in s.ps:
      x=p.exponent
      if x>1: r+='%d,'%x
      else: ones+=1
    r+='1^%d'%ones
    return r+']'

def get_index_biggest(x):
  mx,imx=-1.0,-1
  for i,z in enumerate(x):
    if z>mx: mx,imx=z,i
  return imx

def Lagarias_inequality_rhs(logn):
  # H(n)+exp(H(n))*log(H(n))
  gamma=0.57721566490153286060651209008240243104
  expgamma=1.78107241799019798523650410310717954916964521430343020
  Hn=gamma+logn # +1/(2n)-...
  logHn=log(Hn)
  expHn=exp(Hn)
  return Hn+expHn*logHn

def Lagarias_inequality_rhs_on_n(logn):
  # H(n)/n+exp(H(n))/n*log(H(n)) = x+y*z
  gamma=0.57721566490153286060651209008240243104
  expgamma=1.78107241799019798523650410310717954916964521430343020
  n=exp(logn) # FIXME avoid this!
  x=(expgamma+logn+0.5/n-1.0/(12*n*(n+1)))/n # ~ H(n)/n
  y=expgamma*exp(0.5/n) # ~ exp(H(n))/n = exp(gamma+log(n)+1/(2n)-...)/n
  z=log(gamma+logn) # +O(log(1/n))
  print>>stderr,x,y*z
  return x+y*z

def Lagarias_inequality_rhs_on_n(logn): # rougher version
  # H(n)/n+exp(H(n))/n*log(H(n)) = x+y*z
  gamma=0.57721566490153286060651209008240243104
  expgamma=1.78107241799019798523650410310717954916964521430343020
  x=0.0 # ~ H(n)/n
  y=expgamma # ~ exp(H(n))/n = exp(gamma+log(n)+1/(2n)-...)/n
  z=log(gamma+logn) # +O(log(1/n))
  return x+y*z

m,CA=0,2L
ps=Ps()
logn=log(2.0)
sigma=3 # sigma(2)=3, sigma([1^2])=sigma(2^1*3^1)=12, sigma([2,1^1])=sigma(2^2*3^1*5^1)=168, sigma([4,2,1^4])=sigma(2^4*3^2*5^1*7^1*11^1*13^1)=3249792
rho=1.5 # FIXME rational arithmetic!
f=open('b004490.txt','r')
pg1=Primegen()
pg1.next() # 2 already used
pg1.next() # 3 already used

while len(ps)<len(primes):
  print '%3d\tlogn=%.15f\t%s'%(m+1,logn,ps.show_exponent_list(),),
  if sigma<10000000: print '\tsigma=%d'%(sigma,),
  else:              print '\tsigma=big',
  print '\t%.6f'%rho,
  changeable_places=ps.get_changeable_places()
  e=[]
  for i in changeable_places: e.append(ps[i].epsilon)
  j=get_index_biggest(e)
  k=changeable_places[j]
  p=ps[k].p # primes[k]
  logp=log(p)
  logn+=logp
  if k==len(ps)-1: # A - new prime
    ps.update(k,1,log1p(1.0/p/(p+1.0))/logp,p+1)
    ps.ones+=1
    sigma*=p+1 # (p**2-1.0)/(p-1.0)
    #rho*=(1.0-p**(-2))/(1.0-p**(-1)) # *=(p-1/p^a)/(p-1)
    rho*=(p**2-1.0)/float(p**2-p) # *=(p-1/p^a)/(p-1)
    # add a new zero...
    q=pg1.next()
    logq=log(q)
    ps.new_q(P(q,0,log1p(1.0/q)/logq,1.0))
  else: # B, C - old prime boosted
    exponent=ps[k].exponent
    if exponent==1: # B
      exponent=2
      #sigma*=(p**3-1.0)/(p**2-1.0)
      sigma*=(p**2+p+1); sigma/=(p+1)
      #rho*=(1.0-p**(-3))/(1.0-p**(-2))
      rho*=(p**3-1.0)/float(p**3-p)
      ps.ones-=1
      power_sum=ps[k].power_sum+p**exponent
      epsilon=log1p(1.0/power_sum/p)/logp
      ps.update(k,exponent,epsilon,power_sum)
    else: # C
      exponent+=1
      #sigma*=(p**(exponent+1)-1.0)/(p**(exponent)-1.0)
      sigma*=p**(exponent+1)-1
      sigma/=p**(exponent)-1
      #rho*=(1.0-p**(-exponent-1))/(1.0-p**(-exponent))
      rho*=(p**(exponent+1)-1.0)/float(p**(exponent+1)-p)
      power_sum=ps[k].power_sum+p**exponent
      epsilon=log1p(1.0/power_sum/p)/logp
      ps.update(k,exponent,epsilon,power_sum)
  try:
    b004490=int(f.readline().split()[1])
    assert CA==b004490
  except: # file exhausted
    pass
  print '\t%g'%((Lagarias_inequality_rhs_on_n(logn)-rho),)
  CA*=p
  m+=1
  if m>10000: break

f.close()

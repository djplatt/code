# finds a value x_0 such that, for every x>=x_0,
# the number of odd square-free integers <= x
# is (4/(pi*pi)) x + O^*((9/70) sqrt(x))


import numpy as np
def sqfsieve(N):
          sqf=np.ones(N+1,dtype=bool)
          d=2; sqd=d*d;
          while sqd <= N:
            m = sqd;
            while m <= N:
              sqf[m] = False;
              m += sqd;
            d+=1; sqd= d*d;
          return sqf
qrantop = 4342058
c1 = RBF(4/pi**2); c0 = RBF(9/70); c0sq = c0*c0
Q2=0; qranmax=0; sqf = sqfsieve(qrantop)
for n in [1,3..qrantop]:
            if sqf[n]:
              c0sqn = (c0sq*n).lower()
              errQ = Q2-c1*n
              if (errQ*errQ).upper()>c0sqn:
                qranmax=n;
              Q2+=1
              errQ = Q2-c1*n
              if (errQ*errQ).upper()>c0sqn:
                qranmax=n;
qranmax += 1
print(qranmax)
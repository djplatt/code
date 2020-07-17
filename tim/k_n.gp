/* some random multiplicative function that Tim Trudgian needed */


/* k at a prime power */
/* fordiv(a,b,c) does c for b ranging over the divisors of a */
kp(p,pow)={res=0;fordiv(pow,bi,res=res+p^bi);res}

/* k is mult */
k(n)={f=factorint(n); prod(i=1,matsize(f)[1],kp(f[i,1],f[i,2]))}
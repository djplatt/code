/* generate the coefficients for the product of two L-functions
   of quadratic characters, one even, one odd*/

N1=11;
N2=13;

if(kronecker(N1-1,N1)*kronecker(N2-1,N2)!=-1,quit);

M=1024; /* how many coefficients */

foo(n,a,b)={local(res);res=0;fordiv(n,d,res+=kronecker(d,a)*kronecker(n/d,b));res}

coeffs=vector(M,n,foo(n,11,13));

g(u)=2*exp(u/2-2*Pi*exp(u));


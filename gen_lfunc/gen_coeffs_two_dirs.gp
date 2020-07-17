/* generate the coefficients for the product of two L-functions
   of quadratic characters, one even, one odd*/

N1=precprime(1000);
N2=nextprime(N1+1);

while(kronecker(N1-1,N1)*kronecker(N2-1,N2)!=-1,N2=nextprime(N2+1));

cond=N1*N2;

B=256;
hi_i=128; /* last G(u_m) provided */

M=floor(sqrt(cond)*exp((hi_i+0.5)*2*Pi/B));

foo(n,a,b)={local(res);res=0;fordiv(n,d,res+=kronecker(d,a)*kronecker(n/d,b));res}

for(n=1,M,print(foo(n,N1,N2)," 0 0 0 0 0 0 0"));

quit;

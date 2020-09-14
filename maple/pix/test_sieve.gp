/* 10^15 */
x=10^15;
lam=5069816557000022/2^71;
xi=2^20


l2=1.0/lam/sqrt(2);
phi1(t)=0.5*erfc(log(t/x)*l2);
phi(t)=if(t<x,1-phi1(t),-phi1(t));

sieve(A,B)=
{
   local (p,res);
   p=nextprime(A);
   res=0;
   print("Sieving from ",A," to ",B);

   while(p<=B,res=res+phi(p);p=nextprime(p+1));

   return(res);

}

print("Last positive should be  ",sieve(x-2*xi,x));
print("First negative should be ",sieve(x,x+2*xi));

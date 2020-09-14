
num_sieves=70560;
X=10^23;
sw2=num_sieves*2^33;
X0=X-sw2+1;
print("sqrt X0=",floor(sqrt(X0)));
pi_sqrt_X0=12431880422;
lam=6224003264759175.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/8.0;

B=X;
A=X-sw2-sw2+1;

print("X=",X);
print("X0=",X0);
print("A=",A);

a=intnum(t=0.5,1.0,exp(lam^2*t^2/2)*X0^t/t);
print("a=G(1)-G(1/2)=",a);

b=-12410224303.288;
/*
real(-I*intnum(t=0,14,exp(lam^2*(1/2+I*t)^2/2)*X0^(1/2+I*t)/(1/2+I*t)));
*/
print("b=G(1/2+14i)-G(1/2)=",b);

/* sum of zeros */
d=-6425510748.253;
c=-21680940.515;

print("c=G(1/2+inf)-G(1/2+14i)=",c);
print("d= N G(1/2+inf)- sum G(rho) - G(1/2+14i)=",d);

/* sum of phi */
e=348852.390;
print("e=sum phi(p)=",e);

phi(t,x,lam)=1/2*erfc(log(t/x)/sqrt(2)/lam);

phi1(t,x,lam)=if(t<x,1-phi(t,x,lam),-phi(t,x,lam));

n=2;
res=0;
while(2^n<=B,x=nextprime(A^(1/n));\
   while(x^n<=B,res=res+1/n*phi1(x^n,X0,lam);\
      x=nextprime(x+1););\
   n=n+1;);
f=res;
print("f=sum 1/n phi(p^n) n>1 =",f);

g=sum(m=3,n,1/m*primepi(X0^(1/m)))*1.0+pi_sqrt_X0/2.0;
print("g=sigma 1/n pi(X0^(1/n)) n>1 =",g);




/* sum of primes in 10^23-X0 */
h=11444711269535;
print("h=pi(10^23)-pi(X0)=",h);



print("published value=",1925320391606803968923);
print("my calc        =",a-b+2*d-3*c+e+f-g+h-log(2));
print("diff=",1925320391606803968923-(a-b+2*d-3*c+e+f-g+h-log(2)));

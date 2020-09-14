\p100
/* returns G(inf)-G(1/2+It) */
GG(t)=real(intnum(s=0.5,-1,exp(lam*lam/2*(s+I*t)^2)*X0^(s+I*t)/(s+I*t)));

num_sieves=70560;
X=10^23;
sw2=num_sieves*2^33;
X0=X-sw2+1;
print("sqrt X0=",floor(sqrt(X0)));
pi_sqrt_X0=12431880422;
lam=6224003264759175.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/8.0;

print("G(T1)=",GG(1.1e10));

B=X;
A=X-sw2-sw2+1;

print("X=",X);
print("X0=",X0);
print("A=",A);

a=intnum(t=0.5,1.0,exp(lam^2*t^2/2)*X0^t/t);
print("a=G(1)-G(1/2)=",a);

/* from G_14i */
b=-1.24102243032887041609664681649686732173e10;
print("b=G(1/2+14i)-G(1/2)=",b);
print("            (check) ",real(I*intnum(t=0,14,exp(lam*lam/2*(1/2+I*t)^2)*(X0^(1/2+I*t))/(1/2+I*t))));

/* sum of zeros from G1.6 and sum_G*/
d=-6425510748.253; /* - sum rho - G(14) */
c=-21680940.515; /* -G(14) */

print("c=G(1/2+inf)-G(1/2+14i)=",c);
print("                (check) ",GG(14));

print("d= N G(1/2+inf)- sum G(rho) - G(1/2+14i)=",d);

/* sum of phi from phi2.0 and sum_phi*/
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

g=sum(m=3,log(X0)/log(2),1/m*primepi(X0^(1/m)))*1.0+pi_sqrt_X0/2.0;
print("g=sigma 1/n pi(X0^(1/n)) n>1 =",g);

pistar=a-b+2*d-3*c-log(2)+e+f;
print("pi*(X0)=",pistar);

pix0=pistar-g;
print("pi(X0)=",pix0);

/* sum of primes in 10^23-X0 */
h=11444711269535;
print("h=pi(10^23)-pi(X0)=",h);

print("published value=",1925320391606803968923);
print("my calc        =",pix0+h);
print("diff=",1925320391606803968923-(a-b+2*d-3*c+e+f-g+h-log(2)));

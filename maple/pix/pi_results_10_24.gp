\p100

X=10^24;
X0=X;
sw2=352800/2*2^34;
print("sqrt X0=",floor(sqrt(X0)));
pi_sqrt_X0=37607912018;
lam=6273445730170391.0/2^84;
print("lambda=",lam);
print("      =2^",log(lam)/log(2))
lam2=lam*lam/2;

/* from read_zeros */
T1=20950046000.0;
NT1=69778732700;


/* returns G(1/2+It) */
GG(t)=-real(intnum(s=0.5,-1,exp(lam2*(s+I*t)^2)*X0^(s+I*t)/(s+I*t)));

print("G(1/2+T1)=",GG(T1));
B=X+sw2;
A=X-sw2+1;

print("X=",X);
print("A=",A);
print("B=",B);

a=intnum(t=0.5,1.0,exp(lam2*t^2)*X0^t/t);
print("a=G(1)-G(1/2)=",a);

bb=real(I*intnum(t=0,14,exp(lam2*(1/2+I*t)^2)*X0^(1/2+I*t)/(1/2+I*t)));
b=-3.6619604318105e10;
print("b=G(1/2+14i)-G(1/2)=",b);
print("check               ",bb);

/* sum of zeros from G1.6 and sum_G*/
d=-1.5842462676606395e8; /* - sum rho - G(14) */
d_high=0.0;
c=-9.883459626999901185364694e8; /* -G(14) */

print("c=G(1/2+inf)-G(1/2+14i)=",c);
print("check                   ",GG(T1)-GG(14));

print("d= N G(1/2+inf)- sum G(rho) - G(1/2+14i)=",d);

/* sum of phi from phi2.0 and sum_phi*/
e=-8.53200681e5;/*03472579*/
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




/* sum of primes in X-X0 */
h=0;

/* sum of G(T1) */

print("G(1/2+T1i)=",GG(T1));
k=GG(T1)*(NT1*2-3);
print("k=(2N-3)G(1/2+T1i)=",k);
P=18435599767349200867866;
Li(x)=-eint1(-log(x));
print("R(x)           =",sum(n=1,100,moebius(n)/n*Li(X0^(1/n))));
print("published value=",P);
print("my calc        =",a-b+2*d-3*c+e+f-g+h-log(2)-k);
print("diff=",P-(a-b+2*d-3*c+e+f-g+h-log(2)-k));
/*
print("my calc        =",a-b+2*d_high-3*c+e+f-g+h-log(2)-k);
print("diff=",P-(a-b+2*d_high-3*c+e+f-g+h-log(2)-k));
*/

/*
dd=(P-(a-b-3*c+e+f-g+h-log(2)-k))/2;
print("sum of zeros (d) should be ",dd);
print("sum of zeros (d) was       ",d);
*/

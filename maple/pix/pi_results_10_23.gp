\p100

/* returns G(inf)-G(1/2+It) (roughly) */
GG(t)=real(intnum(s=0.5,-1,exp(lam*lam/2*(s+I*t)^2)*X^(s+I*t)/(s+I*t)));


X=10^23;
pi_sqrt_X=12431880460; /* from nth prime page */
lam_n=6224003264759175.0 /* numerator */
lam_d=83;               /* power of 2 for denominator */
lam=lam_n/2^lam_d;
print("lambda=",lam_n,"/2^",lam_d);
print("      =2^",log(lam)/log(2));
print("      =",lam);
num_sieves=70560;

/* from read_zeros */
T1=1.1155646e10;
NT1=36037434430;


sw=num_sieves*2^34;
B=X+sw/2;
A=X-sw/2+1;

print("X=",X);
print("A=",A);
print("B=",B);

/* numerical integration seems ok
   check with pi_integral.gp if worried */

a=intnum(t=0.5,1.0,exp(lam^2*t^2/2)*X^t/t);
print("a=G(1)-G(1/2)=",a);

/* from G_14i */
b=-1.2410224303293e10;
print("b=G(1/2+14i)-G(1/2)=",b);
print("            (check) ",real(I*intnum(t=0,14,exp(lam*lam/2*(1/2+I*t)^2)*(X^(1/2+I*t))/(1/2+I*t))));

/* sum of zeros from G1.6 and sum_G*/
d=-5.382878886e8; /* - sum rho - G(14) */
c=-2.1680976701924217e7; /* -G(14) */

print("c=G(1/2+T1i)-G(1/2+14i)=",c);
print("                (check) ",GG(14));

print("d= N G(1/2+T1i)- sum G(rho) - G(1/2+14i)=",d);

/* sum of phi from phi2.0 and sum_phi*/
e=-87064.24;
print("e=sum phi(p)=",e);

phi(t,x,lam)=1/2*erfc(log(t/x)/sqrt(2)/lam);

phi1(t,x,lam)=if(t<x,1-phi(t,x,lam),-phi(t,x,lam));

n=2;
res=0;
while(2^n<=B,x=nextprime(A^(1/n));\
   while(x^n<=B,res=res+1/n*phi1(x^n,X,lam);\
      x=nextprime(x+1););\
   n=n+1;);
f=res;
print("f=sum 1/n phi(p^n) n>1 =",f);

/* primepi struggles above X^(1/3) */
g=sum(m=3,log(X)/log(2),1/m*primepi(X^(1/m)))*1.0+pi_sqrt_X/2.0;
print("g=sigma 1/n pi(X^(1/n)) n>1 =",g);

k=GG(T1)*(NT1*2-3);
print("k=(2N-3)G(1/2+T1i)=",k);

pistar=a-b+2*d-3*c-log(2)+e+f-k;
print("pi*(X)=",pistar);

pix=pistar-g;
print("pi(X)= ",pix);

print("published value=",1925320391606803968923);
print("my calc        =",pix);
print("diff=",1925320391606803968923-pix);

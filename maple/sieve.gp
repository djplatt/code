/* 10^13 */
X=1e13;
lam=4646592975839.0/2^60;
sw=566986042.0;


/* 10^15 */
X=1e15;
lam=5069816557000022.0/2^71;
sw=2^35;

/* 10^22 */
X=1e22-273747*2^32+1;
lam=7699934548453755.0/2^79;
sw=273747*2^32;

l2=1.0/lam/sqrt(2);
B=X+sw/2;
A=X-sw/2;

phi(t)=0.5*erfc(log(t/X)*l2);

phi1(t)=if(t<X,1-phi(t),-phi(t));

n=2;
res=0;

print("sieving from ",A," to ",B," with lambda=",lam); 
while(2^n<=B,x=nextprime(A^(1/n));\
   while(x^n<=B,res=res+1/n*phi1(x^n);\
      x=nextprime(x+1););\
   n=n+1;);
print("sigma n=2..",n," 1/n phi(p^n)=",res);

x=nextprime(A);
res=0;
count=0;
while(x<=B,res=res+phi(x);count=count+1;if((count%1000000)==0,print(x));x=nextprime(x+1));
print("sigma phi(p)=",res);

base=2^64;
my_write(n,x)=print(n," ",floor(x/base)," ",x%base);

X=10.0^14;
lam=4695320538125.0/2^60;
sw=100000000;
XI=2^32;

B=X+sw/2;
A=X-sw/2;

phi(t,x,lam)=1/2*erfc(log(t/x)/sqrt(2)/lam);

phi1(t,x,lam)=if(t<x,1-phi(t,x,lam),-phi(t,x,lam));

n=1;
res=0;
count=0;
sum1=0;
sum2=0;
t0=A+XI/2;
t1=t0+XI/2;
print("sieving from ",A," to ",X," with lambda=",lam); 
while(2^n<=X,x=nextprime(A^(1/n));\
   while(x^n<=X,\
      if(x>t1,\
         print("count=",count," sum=",sum1," sum2=",sum2);\
         count=0;\
         sum1=0;\
         sum2=0;\
         t0=t0+XI;\
         t1=t1+XI;);\
      count=count+1;\
      del=t0-x;\
      sum1=sum1+del;\
      sum2=sum2+del*del;\
      x=nextprime(x+1););\
   n=n+1;);
print("finished with n=",n," count=",count," sum1=",sum1," sum2=",sum2);

n=1;
count=0;
sum1=0;
sum2=0;
t0=X+XI/2;
t1=t0+XI/2;
print("sieving from ",A," to ",X," with lambda=",lam); 
while(2^n<=B,x=nextprime(A^(1/n));\
   while(x^n<=B,\
      if(x>t1,\
         print("count=",count," sum=",sum1," sum2=",sum2);\
         count=0;\
         sum1=0;\
         sum2=0;\
         t0=t0+XI;\
         t1=t1+XI;);\
      count=count+1;\
      del=t0-x;\
      sum1=sum1+del;\
      sum2=sum2+del*del;\
      x=nextprime(x+1););\
   n=n+1;);
print("finished with n=",n," count=",count," sum1=",sum1," sum2=",sum2);

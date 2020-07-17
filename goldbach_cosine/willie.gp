P6N(P,N)=P/gcd(P,6*N);

W(P,N)=sum(n=1,P,if(gcd(N^2-n^2,P)==1,1,0));
T(P,N,x)=sum(n=1,P,if(gcd(N^2-n^2,P)==1,frac((x-n)/P)-1/2,0));

W1(P,N)={Np=gcd(N,P);p6n=P6N(P,N);print("Np=",Np," P6N=",p6n);res=1;forprime(p=2,Np,if(Np%p==0,res*=p-1));forprime(p=2,p6n,if(p6n%p==0,res*=p-2));return(res)};

Wd(d,P,N)={p6n=P6N(P,N);sum(n=1,P,if((gcd(n,p6n)==d)&&(gcd(N^2-n^2,P)==1),print(n," ",N^2-n^2);1,0))};

Wd1(d,P,N)={res=1;fordiv(gcd(N,P),p,if(isprime(p),res*=p-1));fordiv(P6N(P,N*d),p,if(isprime(p),res*=p-3));return(res);}

S(P,N,x)={p6n=P6N(P,N);sum(n=1,x,if(gcd(N^2-n^2,P)==1,print(n," ",N^2-n^2," ",gcd(n,p6n));1,0))};

Sd(d,P,N,x)={p6n=P6N(P,N);sum(n=1,x,if((gcd(p6n,n)==d)&&(gcd(N^2-n^2,P)==1),print(n," ",N^2-n^2);1,0))};

/* test with theta(N)=0 */
test_it(N,y,x)={P=1;forprime(p=2,sqrt(2*N),P*=p);p6n=P6N(P,N);fordiv(p6n,p,if(isprime(p),if(Sd(p,P,N,y+x)-Sd(p,P,N,y-x)-3*x/P*Wd1(p,P,N)>0,print("failure at p=",p," x=",x," y=",y))));}

/* test with theta(N)=sqrt(N)/log(N)^3 */
test_it1(N,y,x)={P=1;forprime(p=2,sqrt(2*N),P*=p);p6n=P6N(P,N);fordiv(p6n,p,if(isprime(p),if(Sd(p,P,N,y+x)-Sd(p,P,N,y-x)-3*x/P*Wd1(p,P,N)-sqrt(N)/log(N)^3>0,print("failure at p=",p," x=",x," y=",y))));}


test_it1(400,201,200);


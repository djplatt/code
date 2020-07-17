MAX_X=1000;
Th=vector(MAX_X*2);
forprime(p=2,MAX_X*2,lp=log(p);pn=p;while(pn<=MAX_X*2,Th[pn]=lp;pn*=p));
for(i=3,MAX_X*2,Th[i]=Th[i]+Th[i-1]);

II(X)={7*X/3+sum(a=X,X+X-1,Th[a]^2-(a+a+1)*Th[a])/X^2};

II1(X)=sum(a=X,X+X-1,intnum(t=a,a+1,abs(Th[a]-t)^2))/X^2;

II2(X)=sum(a=X,X+X-1,intnum(t=a,a+1,Th[a]^2-2*Th[a]*t+t^2))/X^2;
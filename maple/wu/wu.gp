NP(N,P)=gcd(N,P);

Pn(P,n)=P/gcd(P,n);

WP(P,N,n)=if(gcd((N-n)*(N+n),P)==1,1,0);

prodp1(N)={local(f);f=factor(N);prod(n=1,length(f[,1]),f[n,1]-1)};
prodp2(N)={local(f);f=factor(N);prod(n=1,length(f[,1]),f[n,1]-2)};

test_WP(P,N)={local(lhs,rhs);lhs=sum(n=1,P,WP(P,N,n));rhs=prodp1(NP(N,P))*prodp2(Pn(P,6*N));return(lhs==rhs)};

CP(P,N,k)=sum(n=1,P,if(WP(P,N,n),cos(2*n*k*Pi/P)));

check_SP(N)={local(P,N2);N2=N*N;P=1;forprime(p=2,sqrt(2*N),P*=p);for(n=1,N-2,if(gcd(N2-n*n,P)==1,return(1)));return(0)};


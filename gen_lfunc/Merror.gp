myexp(x)=if(x<-100000,exp(-100000),exp(x));

calc_C(r,alpha)={local(res);res=1;forprime(p=2,r^(1/alpha),res*=prod(j=1,(r-1)/(p^alpha-1),(r+j-1)/(j*p^alpha)));return(res);}

nu_j(j,mus)=if(j<=2,(mus[j]-1/2)/2,(mus[2]-1)/2);
X_u(u,r)=Pi*r*myexp(2*u/r);

Merror(x,M,N,mus,alpha)={local(r,C,nu,c,X,u);u=x-1/2*log(N);r=length(mus);C=calc_C(r,alpha);nu=1/2+2/r*sum(j=1,r,nu_j(j,mus));c=nu+1/2+alpha;X=X_u(u,r);return(C*2^(r/2)*M^(c-1)*myexp(nu*u-X*M^(2/r))*(1+r*M/(2*X*M^(2/r)+2-c*r))*prod(j=1,r,(1+r*nu_j(j,mus)/(X*M^(2/r)))^nu_j(j,mus)));}

bigpi(x)=x/log(x)*(1+1/2762/log(x));


for(r=2,9,mus=vector(r,x,1);for(n=10,41,print(r," ",2^n," ",ceil(bigpi(solve(M=2^(n/2),2^(2*n),Merror(0,M,2^n,mus,0.6)-1e-60))))));

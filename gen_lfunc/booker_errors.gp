r=2; /* degree */
mus=[1.0,2.0];
alpha=1/r;
nuj(j)=if(j<=2,(mus[j]-0.5)/2,(mus[j]-1.0)/2);
nu=0.5+2/r*sum(j=1,r,nuj(j));
N=60; /* conductor */
A=16;
B=256;

Merror(x,M)={local(u,X,c,C,XMr2);u=x-0.5*log(N);X=Pi*r*exp(2*u/r);XMr2=X*M^(r/2);c=nu+0.5+alpha;C=sqrt(3);if(XMr2<=max(c*r/2-1,r/2),print("Failed XM^(r/2) test.");return(0));print("u=",u);print("X=",X);print("c=",c);return(C*2^(r/2)*M^(c-1)*exp(nu*u-X*M^(2/r))*(1+r*M/(2*X*M^(2/r)+2-c*r))*prod(j=1,r,(1+r*nuj(j)/(X*M^(2/r)))^nuj(j)))}

Fhat(x)={local(u,X,c,C);u=x-0.5*log(N);X=Pi*r*exp(2*u/r);c=nu+0.5+alpha;C=sqrt(3);return(zeta((2*X/r)^r)*2^(r/2)*exp(nu*u-X)/(1-exp(2*Pi*A*(1/2-2*X/r)))*prod(j=1,r,(1+r*nuj(j)/X)^nuj(j)))}

Mworst_case(n)=Merror(n*2*Pi/B,1);

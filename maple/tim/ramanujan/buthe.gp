nu(c,al,K)={local(h);h=(1-al)/K;h^2*c*sum(k=1,K,sum(j=1,k,besseli(0,c*sqrt(2*j*h-j^2*h^2))/(2*sinh(c))));}

nul(c,al,K)={local(h);h=(1-al)/K;h^2*c*sum(k=0,K-1,sum(j=0,k,besseli(0,c*sqrt(2*j*h-j^2*h^2))/(2*sinh(c))));}

mu(c,al,K)={local(h);h=(1-al)/K;h*c*sum(k=1,K,besseli(0,c*sqrt(2*k*h-k^2*h^2))/(2*sinh(c)));}

e1(c,al,K,T,x)={local(eps,x0,vc,uc,B0);eps=c/T;print("eps set to ",eps);x0=x*exp(-al*eps);print("log x0 set to ",log(x0));vc=nu(c,al,K);uc=mu(c,al,K);B0=eps*x0*vc/(2*uc+al);exp(2*eps)*log(exp(eps)*x0)*((2*eps*vc)/log(B0)+2.01*eps/sqrt(x0)+log(log(2*x0^2))/(2*x0));}

e2(c,al,K,T,x)={local(eps,x0,vc,uc,B0);eps=c/T;x0=x*exp(-al*eps);0.16*(1+1/x0)/sinh(c)*exp(0.71*sqrt(c*eps))*log(c/eps);}
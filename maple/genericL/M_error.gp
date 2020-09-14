r=6;
C=22.05;
alpha=1;
N=9449;
mus=vector(r,d,if(d<4,0.5,1.5));
nus=vector(r,n,if(n<=2,(mus[n]-1/2)/2,(mus[n]-1)/2));
nu=1/2+2/r*sum(j=1,r,nus[j]);

k0=r*(nu+alpha+1/2)/2;
k1(x)=Pi*r*exp(2*x/r);
k2=2^((r-2)/2)*r*C*N^((2*alpha+1)/4);

M_error(x,M)=k2*exp(nu*x)*k1(x)^(-k0)*incgam(k0,k1(x)*(M/sqrt(N))^(2/r))*prod(j=1,r,(1+nus[j]*r/k1(x)*(sqrt(N)/M)^(2/r))^nus[j]);

gammar(s)=Pi^(-s/2)*gamma(s/2);
g(u,t)=1/(2*Pi)*exp(-u*I*t)*prod(j=1,r,gammar(1/2+I*t+mus[j]));
G(u,T)=2*intnum(t=0,T,real(g(u,t)));
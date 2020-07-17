/* produce G files for Odd Dirichlet Character L-function */

\p 200

N=11;
r=1;
mus=[1];
m=0;
C=1;
alpha=0;
et=0;
delta=Pi/2*(1-et);
nus=vector(r,n,(real(mus[n])-1)/2);
mu=-1/2+1/r*(1+sum(j=1,r,mus[j]));
K=2*sqrt(2^(r+1)/r*exp(delta*(r-1))/delta)*exp(-Pi*r*et*imag(mu)/4);
X(x)=Pi*r*delta*exp(-delta)*(exp(x)/sqrt(N))^(2/r);
c=real(mu)+1/2+alpha;
c_prime=max(c*r/2-1,0);
lem54(M,x)={local(XX);XX=X(x);if(XX*M^(2/r)<=max(c_prime,r),print("M too small.");0,K*r/2*(exp(x)/sqrt(N))^real(mu)*C*M^c*exp(-XX*M^(2/r))/(XX*M^(2/r)-c_prime)*prod(j=1,r,(1+r*nus[j]/(XX*M^(2/r)))^nus[j]))};

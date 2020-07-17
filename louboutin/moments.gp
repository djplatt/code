/* compute one value */
g(n,p)=n*exp(-n^2*Pi/p);

/* compute chi(n) for the m'th character mod p with prim root pr */
chi(n,m,p,pr)=if(gcd(n,p)>1,0,exp(2*Pi*I*znlog(n,Mod(pr,p))*m/(p-1)));

/* stop when terms decay to exp(-target) */
target=100;

/* compute theta(1,chi) for odd chi as a vector */
moments(p)={local(v,len,pr,ptr);pr=znprimroot(p);if((p%4)==3,len=(p+1)/4,len=(p-1)/4);v=vector(len);ptr=1;for(c=1,len,v[ptr]=sum(j=1,sqrt(target*p/Pi),chi(j,2*c-1,p,pr)*g(j,p));ptr++);return(v)};

/* print the moments */
print_moments(p,max_k=9)={local(extra,len,V);V=moments(p);extra=if((p%4)==3,1,2);len=length(V);for(k=1,max_k,printf("%lu: %10.8e\n",k,2*sum(n=1,len-1,abs(V[n])^(2*k))+extra*abs(V[len])^(2*k)))};


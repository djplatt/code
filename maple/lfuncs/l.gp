gamr(s)=gamma(s/2)*Pi^(-s/2);
Lam(s,L,N,eps,mus)=lfun(L,s)*N^((s-1/2)/2)*eps*prod(n=1,length(mus),gamr(s+mus[n]));
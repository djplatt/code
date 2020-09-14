read("rhos.gp")

lambda=7699934548453755.0/2^79;

x0=10^22-273747*2^32+1;
logx0=log(x0);
e1=sqrt(x0)*exp(lambda^2/8);
l2=lambda^2/2;

G(t1,t2)=e1*real(I*intnum(t=t1,t2,exp(l2*(I*t-t^2)+I*t*logx0)/(0.5+I*t)));

print(G(rhos[1][1],rhos[2][1]));
print(2*G(rhos[2][1],rhos[3][1]));


print(4519*G(rhos[4519][1],rhos[4520][1]));
print(4520*G(rhos[4520][1],rhos[4521][1]));

res1=0;
res2=0;
for(n=1,4521,{res=G(rhos[n][1],rhos[n+1][1]);res1=res1+n*res;res2=res2+res;});
print("sum G(rho)=",res1);
print("G(end)-G(14)=",res2);

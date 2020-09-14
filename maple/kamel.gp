k2(r)=1/sin(Pi/(2*(r+1)));
k3(r)=sqrt(3)/(2*sin(Pi/(3*(r+1))));

k23=k2(3);
k33=k3(3);

rho3(u)=u^(k33-1)*if(u<=2,(1/gamma(k33)-k33*intnum(t=1,u,1/t*(1-1/t)^(k33-1))),if(u<=3,rho3(2)-k33*intnum(t=2,u,rho3(t-1)/t^k33),if(u<=4,rho3(3)-k33*intnum(t=3,u,rho3(t-1)/t^k33),rho3(4)-k33*intnum(t=4,u,rho3(t-1)/t^k33))));

f(t)=1/t*(1-1/t)^(k3(3)-1);

for(t=10,20,print(t/10.0,",",f(t/10.0)));

quit;


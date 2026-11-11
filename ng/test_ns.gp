\p200

H=2089.0/16384.0;
A=4096/21;
t0=10001.7;
z=t0-I/2;



Lam(t)=(1/2-I*t)*Pi^(-I*t/2)*gamma(1/4+I*t/2)*zeta(1/2+I*t);
Lam1(t)=(1/2-I*t)*Pi^(-I*t/2)*gamma(1/4+I*t/2);

W(t)=Lam(t)*exp(Pi*t/4-Pi*(t-t0)^2/H^2);
W1(t)=Lam1(t)*exp(Pi*t/4-Pi*(t-t0)^2/H^2);

k0=round(t0*A);

mysinc(z)=sinc(Pi*z);

this_W=sum(k=-350,350,kk=k0+k;W(kk/A)*mysinc(A*z-kk));

print(this_W/W1(z));

print(zeta(1+I*t0));





\p200

om=727.951335;
om=727.951332982;

A=30610046000;
al=2e18;
et=2*A/al;
li(x)=-real(eint1(-log(x)));
K(y)=sqrt(al/2/Pi)*exp(-al*y^2/2);

H1()=-1/2*intnum(u=om-et,om+et,K(u-om)*u*li(exp(u/2))/exp(u/2));

print(H1());

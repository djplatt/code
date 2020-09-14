print("Tried to use Section 4 of Pintz paper for Ramanujan. Failed.");

e0_1000=1.2e-5;
e0_2000=8.35e-10;
e0_3000=4.51e-13;
li(t)=real(-eint1(-log(t)));
E(x,lx)=7.6+x*e0_3000/lx+e0_3000*(li(x)-1/lx-li(exp(3000)+1/3000))+e0_2000*(li(exp(3000))-1/3000-li(exp(2000)+1/2000))+e0_1000*(li(exp(2000))-1/2000-li(exp(1000)+1/1000));

pil(x,lx)=li(x)-E(x,lx);
pih(x,lx)=li(x)+E(x,lx);

ram(lx)={local(x);x=exp(lx);pih(x,lx)^2-exp(1)*x/lx*pil(x/exp(1),lx-1)};
myli(x)=real(-eint1(-x))/exp(x);
om=727.95;
et=1e-6;
al=1e16;
K(y)=sqrt(al/2/Pi)*exp(-y^2*al/2);

res=intnum(u=om-et,om+et,u*myli(u/2)*K(u-om));

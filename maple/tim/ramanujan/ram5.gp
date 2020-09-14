Li(x)=real(-eint1(-log(x)));
pip(x,k)=Li(x)+sqrt(x)*log(x)*k;
pim(x,k)=Li(x)-sqrt(x)*log(x)*k;
g(x,k)=pip(x,k)^2-exp(1)*x/log(x)*pim(x/exp(1),k);
e0(x,R)={local(X);X=sqrt(log(x)/R);return(sqrt(8*X/Pi)*exp(-X))}

pip1(x,k)=Li(x)+k*x/log(x);
pim1(x,k)=Li(x)-k*x/log(x);

gg(x,k)=pip1(x,k)^2-exp(1)*x/log(x)*pim1(x/exp(1),k);

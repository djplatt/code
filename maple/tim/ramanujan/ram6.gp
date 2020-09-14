g(x)=pip(x)^2-exp(1)*x/log(x)*pim(x/exp(1));

li(x)=-real(eint1(-log(x)));

R0=5.573412;

eps1(lx)=sqrt(8/(Pi*R0))*lx^(1/4)*exp(-sqrt(lx/R0));

/* interpolate Faber & Kadiri's table for |psi(x)-x|/x between exp(50) and exp(1000) (say) */ 
eps(lx)=1.6351e-9-(lx-100)/100*0.047e-9;

therr(x)=eps(log(x))*x;

psith(x)=(1+7.5e-7)*sqrt(x)+3*x^(1/3);

C0=intnum(t=2,exp(50),eps(log(t))/log(t)^2);

pip(x)=li(x)+therr(x)/log(x)+2/log(2)+2*C0+intnum(t=exp(50),x,eps(log(t))/log(t)^2);

pim(x)=li(x)-(therr(x)+psith(x))/log(x)+2/log(2)-2*C0-intnum(t=exp(50),x,(t*eps(log(t))+psith(t))/(t*log(t)^2));

duspim(x)=x/log(x)*(1+1/log(x)+2/log(x)^2);
duspip(x)=x/log(x)*(1+1/log(x)+2.334/log(x)^2);


pix_xe(x)={local(lx,lxe,xe);lx=log(x);lxe=lx-1;xe=x/exp(1);return(li(x)-li(xe)+eps(lx)*x/lx+eps(lxe)*xe/lxe+intnum(t=xe,x,eps(log(t))/(log(t)^2)));}

g1(x)=pix_xe(x)^2-pim(x/exp(1))*(exp(1)*x/log(x)-pix_xe(x)-pip(x));

pix_xe1(x)={local(lx,lxe,xe);lx=log(x);lxe=lx-1;xe=x/exp(1);return(x/lx*(1+1/lx+2.334/lx^2)-xe/lxe*(1+1/lxe+2/lxe^2));}

pix_xe2(x)=(1+eps(log(x)))*(li(x)-li(x/exp(1)))+2*eps(log(x))*(x/log(x)+x/(exp(1)*log(x/exp(1))));

g2(x)=pix_xe2(x)^2-pim(x/exp(1))*(exp(1)*x/log(x)-pix_xe2(x)-pip(x));


trudpp(x)=li(x)+1.001*sqrt(8/(Pi*R0))*x*(log(x))^(-0.75)*exp(-sqrt(log(x)/R0));
trudpm(x)=li(x)-1.001*sqrt(8/(Pi*R0))*x*(log(x))^(-0.75)*exp(-sqrt(log(x)/R0));

trudg(x)=trudpp(x)^2-exp(1)*x/log(x)*trudpm(x/exp(1));
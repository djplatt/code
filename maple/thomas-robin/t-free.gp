R(t,n)={local(res,tt);tt=t*1.0;res=1.0;fordiv(n,p,if(isprime(p),res*=(1-p^-tt)/(1-1.0/p)));return(res/log(log(n)))}

pn_lo(n)={local(ln,lln);ln=log(n);lln=log(ln);n*(ln+lln-1+(lln-2.1)/ln)};
pn_hi(n)={local(ln,lln);ln=log(n);lln=log(ln);n*(ln+lln-1+(lln-2)/ln)};

sum1_hi(pn_h)={local(lx);lx=log(pn_h);exp(Euler)*lx*(1+0.006/lx^3)};

sum2_hi(pn_l)=exp(2/pn_l);

log_Nn_lo(n)={local(ln,lln);ln=log(n);lln=log(log(n));n*(ln+lln-1+(lln-2)/ln)};

g(t,n)={local(pn_l,pn_h);pn_l=pn_lo(n);pn_h=pn_hi(n);sum1_hi(pn_h)*sum2_hi(pn_l)/zeta(t)/log(log_Nn_lo(n))/exp(Euler)};

gB(t,pn)={local(lpn,spn);lpn=log(pn);spn=sqrt(pn);exp(2/pn)*lpn/(zeta(t)*(1-(3*lpn+5)/(8*Pi*spn))*log(pn-spn*lpn^2/(8*Pi)))};

ginf(t,pn)={local(lpn);lpn=log(pn);exp(2/pn)*log(pn)*(1+0.2/lpn^3)/(zeta(t)*log(pn-0.5*pn/lpn^3))};

B=solve(x=1e10,1e30,4.92*sqrt(x/log(x))-3e12);

gB1(t,pn)={local(lpn,spn);lpn=log(pn);spn=sqrt(pn);exp(2/pn)*lpn*exp(0.00559/lpn^2)/(zeta(t)*log(pn-spn*lpn^2/(8*Pi)))};

pn=precprime(solve(t=1e10,1e100,t+0.5*t/log(t)^3-log(10^(10^13.1))));

f(t,et,k)=(1.02/(t-1)/log(t)+et/k/log(t)^k+(k+2)*et/(k+1)/log(t)^(k+1));

gB2(t,pn,et,k)={local(lpn,spn);lpn=log(pn);spn=sqrt(pn);exp(2/pn)*lpn*exp(f(pn,et,k))/(zeta(t)*log(pn-spn*lpn^2/(8*Pi)))};

C1=intnum(t=B,+oo,(1.388e-10*t+1.4262*sqrt(t))*(1+log(t))/t^2/log(t)^2)

gB3(t,pn,B)={local(lpn,spn);lpn=log(pn);spn=sqrt(pn);exp(2/pn)*lpn*exp(1.02/(pn-1)/lpn+lpn/spn/8/Pi+C1+((lpn+3)*sqrt(B)-(log(B)+3)*sqrt(pn))/(4*Pi*sqrt(pn*B)))/(zeta(t)*log(pn-spn*lpn^2/(8*Pi)))};

ginf1(t,pn)={local(lpn);lpn=log(pn);exp(2/pn)*log(pn)*exp(f(pn,0.5,3))/(zeta(t)*log(pn-1.4262*sqrt(pn)-1.338e-10*pn))};

pi_low(x)=if(x<1000000,primepi(x),x/log(x)*(1+1/log(x)+1.8/(log(x))^2));

pi_high(x)=if(x<1000000,primepi(x),x/log(x)*(1+1/log(x)+2.51/(log(x))^2));

phi(t,lam,x0)=1.0-0.5*erfc(log(t/x0)/lam/sqrt(2));

myexp(t)=if(t<-100000,exp(-100000),exp(t));

phi1(t,lam,x0)={
/*print("t=",t);*/
-myexp(-log(t/x0)^2/2/lam^2)/lam/t/sqrt(2*Pi);
}

lam=6273445730170391*(2^-84);

x0=10^24;

w=6061057848115200/2;

phi_err(x,lam,x0)=phi(x,lam,x0)*pi_high(x)+intnum(t=1,x,phi1(t,lam,x0)*pi_low(t))

phi_m_err(m,x,lam,x0)=intnum(t=1,x^(1/m),t^(m-1)*phi1(t^m,lam,x0)*pi_low(t^m));

all_phi_err(w,lam,x0)={
   local (z,x);
   x=x0-w;
   z=phi(x,lam,x0)*pi_high(x);

   return(sum(m=1,floor(log(x0-w)/log(2)),z/m+phi_m_err(m,x,lam,x0)));
}

phi_up(t,lam,x0)=0.5*erfc(log(t/x0)/lam/sqrt(2));

myexp_up(t)=if(t<-100000,0,exp(t));

phi1_up(t,lam,x0)={
/*print("t=",t);*/
myexp_up(-log(t/x0)^2/2/lam^2)/lam/t/sqrt(2*Pi);
}


phi_m_up_err(m,x,lam,x0)=intnum(t=x^(1/m),[1],t^(m-1)*phi1_up(t^m,lam,x0)*pi_high(t^m));

all_phi_up_err(M,w,lam,x0)={
   local (z,x);
   x=x0+w;
   z=phi_up(x,lam,x0)*pi_low(x);

   return(sum(m=1,M,z/m-phi_m_up_err(m,x,lam,x0)));
}
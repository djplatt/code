
N_adj(t)=
{
   local (lnt);

   lnt=log(t);
   return(0.137*lnt+0.443*log(lnt)+2.462);
}

N_main(t)=
{
   local (t_pi);

   t_pi=t/(Pi*2.0);
   return(t_pi*(log(t_pi)-1.0));
}

N_low(t)=
{
   return(N_main(t)-N_adj(t));
}

N_high(t)=
{
   return(N_main(t)+N_adj(t));
}

myeint1(t)=
{
   if(t>10000,return(1.2e-4347),return(eint1(t)));
}

/* bound for G(sig+I*t) */
G(x,lam,t,sig)=
{
   local (lt2,lam2);

   lam2=lam*lam;
   lt2=(lam2*t*t)/2;
   return(x^sig/(2*lam2)*exp(lam2*sig^2/2)*myeint1(lt2));
}

erf(t)=1-erfc(t);


/* error from taking G(T1)=0 */
err1(x,lam,T1,NT1)=
{
   return(2*NT1*G(x,lam,T1,1/2));
}

/* error from zeros missed above T2 (Re(rho) in (0,1)) */
err2(x,lam,T2,NT2,alp)=
{
   return(2*NT2*G(x,lam,T2,1));/*+exp(lam^2/2)*x*2^((alp-1)/2)*lam^(-alp-1)*incgam((alp+1)/2,lam^2*T2^2/2));*/
}

/* error from zeros between T1 and T2 (with Re=1/2) */
err3(x,lam,T1,NT1,T2,NT2,alp)=
{
   local (tm);
   tm=incgam((alp+1)/2,lam^2*T1^2/2);/*+incgam((alp+1)/2,lam^2*T2^2/2);*/
   tm=tm*2^((alp+1)/2)*lam^(-alp-1);
   return(2*(NT1*G(x,lam,T1,1/2)+NT2*G(x,lam,T2,1/2))+x*exp(lam^2/2)*tm);
}

all_errs(x,lam,T1,NT1,T2,NT2,alp)=err1(x,lam,T1,NT1)+err2(x,lam,T2,NT2,alp)+err3(x,lam,T1,NT1,T2,NT2,alp);

sieve_width(x,lam,err)=
{
   local(t);

   t=lam*(sqrt(2*log(lam*x/err))+1);
   return(exp(t)*x-exp(-t)*x);
}

sieve_lb(x,lam,err)=x*exp(-lam*(sqrt(2*log(lam*x/err))+1));

sieve_ub(x,lam,err)=x*exp(lam*(sqrt(2*log(lam*x/err))+1));

main1()=
{

   logT1=9;

   T1=10^logT1;
   NT1=N_high(T1);
   T2=2.44e12;
   NT2=1e13;
   alp=log(NT1)/log(T1)-1;

   print("Zeros found to height 10^",logT1);
   print("Using about ",floor(NT1)," zeros.");
   print("Taking N(t)<=t^(1+",alp,")");


   x=10^22;
   print("x0=",x);

   count=0;
   while(count<3,
   lam=solve(lam=2^(-50),2^(-10),err1(x,lam,T1,NT1)+err2(x,lam,T2,NT2,alp)-0.124);
   ilam=lam;
   pow=0;
   while(ilam<2^53,ilam=ilam+ilam;pow=pow+1;);
   ilam=ilam/2;
   pow=pow-1;
   ilam=round(ilam);
   print("lambda=",ilam,"*2^-",pow);
   lam=ilam/2^pow;
   print("lambda=2^(",log(lam)/log(2),")");

   dx=1e22-ceil(solve(x=1e21,1e22,sieve_ub(x,lam,1/4)-1e22));
   dx=ceil(dx/2^32);
   print("need to do ",dx," sieves of length 2^32 above x0");
   x=1e22-dx*2^32;
   print("x0 is now ",x);
   count=count+1;);

   dx=x-floor(sieve_lb(x,lam,1/4));
   dx=ceil(dx/2^32);
   print("need to do ",dx," sieves of length 2^32 below x0");

   print("total error from zeros=",all_errs(x,lam,T1,NT1,T2,NT2,alp));
   print("sieve error < 1/4.");

}

main()=
{
   local (x,lam,T1,alp,e1,e3);

   T1=1073546000;
   NT1=N_low(T1);
   T2=2.44e12;
   NT2=1e13;
   alp=log(NT1)/log(T1)-1;

   print("Zeros found to height ",T1);
   print("Using about ",floor(NT1)," zeros.");
   print("Taking N(t)<=t^(1+",alp,")");

   asw=273747*2^32*2;
   x=10^22-asw/2;
   print("x=",x);

/*
   lam=solve(l=2^(-50),2^(-10),err1(x,lam,T1,NT1)-0.249);
   print("lambda=2^(",log(lam)/log(2),")");
*/

   lam=7699934548453755.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/512.0;
/*
   lam=lam/1.081;
*/
   print("lambda=2^(",log(lam)/log(2),")");

   e1=err1(x,lam,T1,NT1);

   print("err1 (taking G(1/2+iT1)=0) =",e1);

   e3=err3(x,lam,T1,NT1,T2,NT2,alp);
   print("err3 (zeros missed between T1 and T2)=",e3);

   e2=err2(x,lam,T2,NT2,alp);
   print("err2 (zeros missed beyond T2) =",e2);

   /*print("error from T1 alone=",err2(x,lam,T1,NT1,alp));*/

   print("total error=",e1+e2+e3);

   serr=solve(serr=0.000001,1.0,sieve_width(x,lam,serr)-asw);

   print("sieving error=",serr);
   print("actual sieve width =",asw);

   print("expected number of primes in sieve=",round(asw/log(x)/1000.)*1000);
   print("maximum y for which pi(y) will be known=",x+asw/2);
}  

/*
lam=7699934548453755.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/512.0;
asw=273747*2^32*2;
x=10^22-asw/2;
lam2=lam*lam/2.0;

phi_hat(s)=exp(lam2*s*s)*x^s/s;

rho1=14.13472514173469379045725198356247027078425711569924317568556746014996342981;
rho2=21.022039638771554992628479593896902777334340524902781754629520403587;
rho3=25.01085758014568876321379099256282181865954967255799667249654200674;

G_14()=intnum(t=0.0,14.0,phi_hat(0.5+I*t));

x=1000000;
lam2=1/1024.;
lam2=lam2*lam2/2.0;
G_A()=intnum(s=0.5,1,phi_hat(s));

print("int from 14 to rho1=",intnum(t=14.0,rho1,phi_hat(0.5+I*t)));
print("2*int from rho1 to rho2=",2*intnum(t=rho1,rho2,phi_hat(0.5+I*t)));
print("3*int from rho2 to rho3=",3*intnum(t=rho2,rho3,phi_hat(0.5+I*t)));
*/

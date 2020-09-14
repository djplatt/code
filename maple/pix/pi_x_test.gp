
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


main()=
{
   local (x,lam,T1,alp,e1,e3);

   T1=2546000;
   NT1=N_low(T1);
   T2=2.44e12;
   NT2=1e13;
   alp=log(NT1)/log(T1)-1;

   print("Zeros found to height ",T1);
   print("Using about ",floor(NT1)," zeros.");
   print("Taking N(t)<=t^(1+",alp,")");

   x=10^13;

   lam=solve(l=2^(-50),2^(-10),err1(x,l,T1,NT1)-0.24);
   print("lambda=",round(lam*2^60),"/2^60");
   lambda=round(lam*2^60)/2^60;
   print("lambda=2^(",log(lam)/log(2),")");

   e1=err1(x,lam,T1,NT1);

   print("err1 (taking G(1/2+iT1)=0) =",e1);

   e3=err3(x,lam,T1,NT1,T2,NT2,alp);
   print("err3 (zeros missed between T1 and T2)=",e3);

   e2=err2(x,lam,T2,NT2,alp);
   print("err2 (zeros missed beyond T2) =",e2);

   /*print("error from T1 alone=",err2(x,lam,T1,NT1,alp));*/

   print("total error=",e1+e2+e3);

   asw=sieve_width(x,lam,0.5);

   print("sieving error=",0.5);
   print("actual sieve width =",asw);

   print("expected number of primes in sieve=",round(asw/log(x)/1000.)*1000);
   print("maximum y for which pi(y) will be known=",x+asw/2);
}  


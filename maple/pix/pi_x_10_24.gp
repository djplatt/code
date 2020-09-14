
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
   return(x^sig/(2)*exp(lam2*sig^2/2)*myeint1(lt2));
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

main1(x,T1)=
{
/*
   T1=1.1e10;
*/
   NT1=N_high(T1);
   T2=2.44e12;
   NT2=1e13;
   alp=log(NT1)/log(T1)-1;

   print("Zeros found to height ",T1);
   print("Using about ",floor(NT1)," zeros.");
   print("Taking N(t)<=t^(1+",alp,")");



   print("x=",x);

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

   x0=floor(solve(w=x/2,x,w+sieve_width(w,lam,0.25)-x));
   if((x0&1)==0,x0=x0+1);
   print("Centring at x0=",x0);
   asw=x-x0;
   asw=ceil(asw/2^34)*2^34;
   print("sieve width =",asw);
   print("We need ",ceil(asw/2^32/2)," 2^32 sieves each side.");

   print("total error from zeros=",all_errs(x0,lam,T1,NT1,T2,NT2,alp));
   sieve_err=solve(err=0.1,0.25,sieve_width(x0,lam,err)-asw);
   print("total sieving error=",sieve_err);

}


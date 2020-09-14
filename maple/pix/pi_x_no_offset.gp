small_pi(n)=if(x<32299,0,n/log(n)*(1+1/log(n)+1.8/(log(n))^2));
big_pi(n)=if(x<355991,0,n/log(n)*(1+1/log(n)+2.51/(log(n))^2));
phir(t,x,lam)=1/2*erfc(log(t/x)/sqrt(2)/lam);

phi(s,x,lam2)=x^s*exp(lam2*s*s)/s;

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
/* version up sig line */ 
G(x,lam,t,sig)=
{
   local (lt2,lam2);

   lam2=lam*lam;
   lt2=(lam2*t*t)/2;
   return(0.5*x^sig*exp(lam2*sig^2/2)*myeint1(lt2));
}

/* version up -1 line 
   This is B(sig,T) in my thesis */

GA(x,lam,t,sig)=
{
  return(exp(lam^2/2*(1-t^2))*(1/x/t^2/lam^2+x^sig/log(x)/t));
} 

B_dash(x,lam,t,sig)=
{
   local (lam2,t2);
   lam2=lam*lam;
   t2=t*t;

   return(-exp(lam2*(1-t2)/2)*(x^sig/log(x)*(lam2*t2+1)/t2+1/lam2/x*(lam2*t2+1)/t2/t));
}

err1B(x,lam,T1,T2)=
{
  return(2*(GA(x,lam,T2,0.5)*N_high(T2)-GA(x,lam,T1,0.5)*N_high(T1)-intnum(t=T1,T2,B_dash(x,lam,t,0.5)*N_high(t))));
}

erf(t)=1-erfc(t);


/* error from taking G(T1)=0 */
err3(x,lam,T1,NT1)=
{
   return(2*NT1*G(x,lam,T1,1/2));
}

err3A(x,lam,T1)=
{
   return(2*N_high(T1)*GA(x,lam,T1,-1.0));
}

myexp(x)=if(x>1000,exp(1000),if(x<-1000,exp(-1000),exp(x)));

/* Error from real part of Phihat for zeros with
   imaginary part > T assuming real part = sigma */

calc_alp(T)=log(N_high(T))/log(T);

G_sig(x,lam,T,sig)=
{
   local (alp,l2);
   l2=lam*lam;
   alp=calc_alp(T);
   return(GA(x,lam,T,sig)*((l2*T*T+2)/(l2*T^(2-alp))-N_low(T)));
}

/* errors from zeros above T2, no RH here */
err2A(x,lam,T2)=G_sig(x,lam,T2,1)+G_sig(x,lam,T2,0);

/* error from zeros above T1, RH here */
err1A(x,lam,T1)=2*G_sig(x,lam,T1,0.5);


/* error from zeros missed above T2 (Re(rho) in (0,1)) */
err2(x,lam,T2,NT2)=
{
   local(alp,lam2,lt2);

   lam2=lam*lam;
   lt2=lam2*T2*T2/2;
   alp=log(T2)+lt2-log(NT2)-log(2);
   if(alp<0,return(2.0));
   return((x+0.5)*exp(lam2/2)*(myeint1(lt2)+alp*myexp(-alp*T2)));
}



/* error from zeros between T1 and T2 (with Re=1/2) */
err1(x,lam,T1,NT1)=
{
   local (alp,lam2,lt2);

   lam2=lam*lam;
   lt2=lam2*T1*T1/2;
   alp=log(T1)+lt2-log(NT1)-log(2);
   if(alp<0,return(2.0));
   return(2*sqrt(x)*exp(lam2/8)*(myeint1(lt2)*NT1+alp*myexp(-alp*T1)));
}

min_one_err(lam,x)=exp(lam^2/2)/(2*Pi*x*lam)*(5*sqrt(2*Pi)+2/lam);

all_errs(x,lam,T1,NT1)=min_one_err(lam,x)+err1(x,lam,T1,NT1)+err2(x,lam,T2,NT2)+err3(x,lam,T1,NT1);

all_errsA(x,lam,T1)=min_one_err(lam,x)+err1A(x,lam,T1)+err2A(x,lam,T2)+err3A(x,lam,T1);

all_errsB(x,lam,T1,T2)=min_one_err(lam,x)+err1B(x,lam,T1,T2)+err2A(x,lam,T2)+err3A(x,lam,T1);

calc_tau(x,lam,err)=lam*(sqrt(2*log(lam*x/err))+1);

sieve_width(x,lam,err)=
{
   local(t);

   t=calc_tau(x,lam,err);
   return(x*(exp(t)-exp(-t)));
}

sieve_lb(x,lam,err)=x*exp(-lam*(sqrt(2*log(lam*x/err))+1));

sieve_ub(x,lam,err)=x*exp(lam*(sqrt(2*log(lam*x/err))+1));

set_lam(x,T1,NT1,err)=
{
   local(t,sx);
   sx=2*sqrt(x);
   t=T1*T1/2;
   return(solve(lam=2^(-100),2^(-10),sx*myeint1(lam*lam*t)*NT1-err));
}

evenp(x)=2.0*(x/2.0)==x;


main1A(x,T1,NT1,T2,NT2,nits)=
{
/*
   T1 is the height to which we have zeros
   nits is number of 2^34 sieves we can do in one batch job
     e.g. 7 cores * 240 sieves in about 10 days at 10^23
   NT1 is the number of zeros we have
   T2 is the height to which RH is known
   NT2 is the number of zeros to height T2
*/

   print("Assuming 1st ",NT2," zeros respect RH.");
   print("  ie RH holds to t=",T2);

   print("Zeros found to height ",T1);
   print("Using about ",floor(NT1)," zeros.");
   /*
   print("Taking N(t)<=t^(1+",alp,")");
   */


   print("x=",x);

   all_error=1.0;
   print("Sharing out total error of ",all_error);

   print("*******************computing lam*********************");
   lam=solve(lam=2^(-100),2^(-10),all_errsB(x,lam,T1,T2)-all_error/2.0);
if(lam*x<=all_error/2.0,print("Error too large................"));
/*
   print("********************fixing lam***********************");
   lam=6273445730170391*2^(-84);
*/
   ilam=lam;
   pow=0;
   while(ilam<2^53,ilam=ilam+ilam;pow=pow+1;);
   ilam=ilam/2;
   pow=pow-1;
   ilam=round(ilam);
   print("lambda=",ilam,"*2^-",pow);
   lam=ilam/2^pow;
   print("lambda=2^(",log(lam)/log(2),")");

   sieve_err=all_error/2.0;
   ta=calc_tau(x,lam,sieve_err);
   asw=x*(exp(ta)-exp(-ta));
print("Computed width="asw);
   x0=x;
/*
   print(asw/2^34);
*/
   num_sieve_its=ceil(asw/2^34/nits/2)*nits*2;
   asw=num_sieve_its*2^34;
   print("Sieve start =",x-asw/2+1);
   print("sieve width =",asw);
   print("We need ",num_sieve_its," 2^34 sieves in total.");
   print("i.e. we need ",asw/2^34/nits," jobs, each doing ",nits," iterations.");
   x0=x;
   print("final x0=",x0);
   print("final sieve start =",x-asw/2+1);
   print("      sieve end   =",x+asw/2);

   print("total error from zeros=",all_errsB(x0,lam,T1,T2));
   print("   zeros error 1=",err1B(x0,lam,T1,T2));
   print("   zeros error 2=",err2A(x0,lam,T2));
   print("   zeros error 3=",err3A(x0,lam,T1));
   print("   -1 int error =",min_one_err(lam,x0));

   sieve_err=solve(err=1e-20,all_error,sieve_width(x0,lam,err)-asw);

   print("total sieving error=",sieve_err);

   print("total error="sieve_err+all_errsB(x0,lam,T1,T2));
/*
   print("main zeros error term=",2*sqrt(x0)/log(x0)*intnum(t=T1,T2,exp(-t^2*lam^2/2)*log(t)/2/Pi/exp(1)));

   I1=intnum(t=0.5,1.0,phi(t,x0,lam*lam/2));

   print("Int from 1/2 to 1 =",I1);
*/

}

/*
   sieve2.0 takes 80 mins at x=10^24 (7 per node)
   so 240 its = 320 hours = 13 1/3 days (10% leeway)
*/

x=10^24; 
/* T1=29988446000; /* Maximum height to which we have zeros */
T1=20950046000;
NT1=69778732700;
T2=30610046000;
NT2=103800788359;

main1A(x,T1,NT1,T2,NT2,240*7);


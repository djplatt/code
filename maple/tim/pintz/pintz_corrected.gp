R=5.573412;

simp(lx)={local(del);del=2/sqrt(R*lx);return((844+R^2)*(3)^(8*del/3)*4^del*(lx/R)^(2+del)*exp(-2*sqrt((1-8*del/3)*lx/R))*(1+0.67*sqrt(R/lx)+R/21212/lx))};

myroundup(x,n)={local(e);e=floor(log(x)/log(10))-n;return(ceil(x/10^(e))*10^(e))};

/*
Error from Lemma 1
*/
t1(lx,T)=2.0*lx^2/T;

/*
Error from (12)
*/
t2(lx,T,del,H)=exp(-lx*del)*(1/2/Pi*log(T/2/Pi)^2+1.8642);

KK=0;llam=0;

/* Number of zeros in box. We were using */

N(sig,T)=17.253*T^(8*(1-sig)/3)*log(T)^(5-2*sig)+3*log(T)^2;

/* but could easily use but use

N(sig,T)=2.177*log(T/2)^(2*sig)*log(T)^(5-4*sig)*T^(8/3*(1-sig))+5.633*log(T)^2;
*/

t3c(lx,T,del,H,K,lam)={return(2*lam/T*sum(k=0,K-1,exp(k*log(lam)-lx/(R*(log(T)-k*log(lam))))*N(1-del,T/lam^k)));};

/*
t3b(lx,T,del,H,K,lam)=2*lam/T*(6*log(T)^2*sum(k=0,K-1,exp(k*log(lam)-lx/(R*(log(T)-k*log(lam)))))+2.2*T^(8*del/3)*log(T)^(3+2*del)*sum(k=0,K-1,exp(k*log(lam)-lx/(R*(log(T)-k*log(lam)))-8*del/3*k*log(lam))));
*/

t3a(lx,T,del,H)={local(this_e,best_e,K,lam);K=1;lam=(T/H);best_e=t3c(lx,T,del,H,K,lam);while(true,K++;lam=(T/H)^(1/K);this_e=t3c(lx,T,del,H,K,lam);if(this_e>best_e,KK=K-1;llam=(T/H)^(1/KK);return(best_e),best_e=this_e))};

err(lx,lx1,T,del,H)={local(tt1,tt2,tt3);tt1=t1(lx1,T);tt2=t2(lx,T,del,H);tt3=t3a(lx,T,del,H);return(tt1+tt2+tt3)};

/*
olderr(lx)=sqrt(8/Pi/R)*lx^0.25*exp(-sqrt(lx/R));

simp_err(lx)=err(lx,exp(2*sqrt(lx/R)),2/sqrt(R*lx),H);
*/

/*
Here we decrease x by 1 since our bounds are only true at half odd integers.
We then iterate over T=10^16..10^43, to find a sweet spot. We fix delta=0.01
so that our bound N(sig,t) works.
*/
chk(llx)={local(lx,lx1,e0,T,lT,del,T0,del0,e1lT0,best_K,best_lam);lx=log(exp(llx)-1);lx1=llx+500;e0=1;T=1e16;lT=16;while(T<1e44,del=0.01;while(del<=0.01,e1=err(lx,lx1,T,del,H);if(e1<e0,e0=e1;T0=T;lT0=lT;del0=del;best_K=KK;best_lam=llam);del+=0.001);T*=10;lT++);printf("$%d$ & $10^{%d}$ & %.3f & %d & %4.3f & %5.4E \\\\\n",lx,lT0,del0,best_K,best_lam,myroundup(e0,4));return(myroundup(e0,4));};

/* RH is true up to */
H=2.5e12;


for(n=0,18,chk(1000+500*n));

/* Gourdon's bound for comparison */
/*
H=2445999556030;

for(n=0,10,chk(1000+500*n));
*/


/*
old_N(sig,T)=2*T*log(1+9.8/2/T*(3*T)^(8*(1-sig)/3)*log(T)^(5-2*sig))+103*log(T)^2;

new_N(sig,T)=9.8*(3*T)^(8*(1-sig)/3)*log(T)^(5-2*sig)+103*log(T)^2;

intver(lx,T,del,H)=2*(exp(-lx/R/log(T))*new_N(1-del,T)/T-exp(-lx/R/log(H))*new_N(1-del,H)/H+intnum(t=H,T,exp(-lx/R/log(t))/t^2*(1-lx/R/log(t)^2)*new_N(1-del,t)));

*/
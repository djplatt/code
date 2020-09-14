/*
   zeta-fft-error.gp

   PARI code to calculate parameters for
   windowed FFT based Zeta calculator
*/

MIN_T0=5000
MIN_T01=1000 /* used in inter error calc because MIN_T0 causes underflow */
/* x is maximum width of approximation = 1/(2B)
   h
   k
*/

lngerr(m,z)=abs(bernfrac(2*m))*2^m/(2*m*(2*m-1)*abs(z)^(2*m-1))

tayerr(x,h,k)=
{
   return(2^((k+5)/2)*Pi^(k+1/2)*h^(k+1)*x^k/gamma((k+2)/2))
}

/* Don't call this with very large sig */

C(sig,t0,h,k)=
{
   local(s3,s2,c1,sm,tm,m);

   s2=(sig-1)/2;
   s3=(sig+0.5)^2/(2*h^2);
   c1=2^((6*k+7-sig)/4)*Pi^(k+0.5)*exp(1/(sig+sig))*gamma(s2+1);
   m=0;
   sm=0.0;
   while(m<=s2,
      tm=(sig+0.5)^k*h^(m+1)*2^(m/2-0.5)*(gamma(m/2+0.5)-incgam(m/2+0.5,s3));
      tm=tm+h^(m+k+1)*2^(m/2+k/2-0.5)*incgam((m+k+1)/2,s3);
      sm=sm+t0^(s2-m)/(gamma(s2-m+1)*gamma(m+1))*tm;
      m=m+1);
   return(c1*sm);
}


/* returns error/k! */
Gtwiderr(A,sig,t0,h,k)=
{
  local(s2,s3);

  s2=sig+sig-1;
  s3=sig+sig+1;

  l=0;
  tm=0.0;
  while(l<=(sig-1)/2,
     tm=tm+(1+1/(A*Pi*(4*l+1)))*((4*l+1)^2+t0^2)^(k/2)/l!*exp((4*l+1)^2/(8*h*h)-A*Pi*(4*l+1)/2);
     l=l+1);
  tm=tm*2^(k+3)*Pi^(k+1);
  tm=exp(log(tm)-MIN_T0*MIN_T0/(2*h^2));

  return((2*(1+1/(A*Pi*s2))*C(sig,t0,h,k)*exp(s3^2/(8*h^2)-A*Pi*s2/2)+tm)/k!);
}


maxGtwiderr(A,t0,h,K)=
{
   local(mx,tp,sig,ttp);
   mx=0;
   for(k=0,K-1,
      sig=3;
      tp=Gtwiderr(A,3,t0,h,k);
      while((sig < 99),
         sig=sig+2;
/*         print("in maxGtwiderr with sig=",sig); */
         ttp=Gtwiderr(A,sig,t0,h,k);
         if(ttp<tp,tp=ttp););
      if(tp>mx,mx=tp));
   return(mx);
}

/* returns error/k! */
gtwiderr(B,h,k)=
{
   local(bh);

   bh=B^2/(8*h^2);

   return(8*(Pi*B)^k*(exp(-bh)+2^(1.5*k-0.5)*(h/B)^(k+1)*incgam(0.5+k/2,bh))/gamma(k+1));
}

maxgtwiderr(B,h,K)=
{
   local(mx,tp,kmax);
   mx=0;
   for(k=0,K-1,
      tp=gtwiderr(B,h,k);
/*      print("g~ error at k=",k," is ",tp); */
      if(tp>mx,
          mx=tp;
/*          kmax=k;    */
          ));
/*   print("worst case k =",kmax); */
   return(mx);
}

/* error taking M terms in sum for f^ */
fhatsumerr(M,sig,t0,h)=
{
   local(s2,resid);

   s2=sig+sig-1;
   resid=2*Pi^1.25*exp(1/(8*h^2)-MIN_T0*MIN_T0/(2*h^2));

   return(resid+C(sig,t0,h,0)*exp(s2^2/(8*h^2))*Pi^(s2*(-0.25))*M^(1-sig)/(sig-1));
}

bestfhatsumerr(M,t0,h)=
{
   local(sig,best_err,err);

   sig=3;
   best_err=fhatsumerr(M,sig,t0,h);
   while(sig<5000,
         sig=sig+100;
         err=fhatsumerr(M,sig,t0,h);
         if(err<best_err,best_err=err););
   return(best_err);
}

fhattwiderr(sig,t0,h,A)=
{
  local(ts);
  ts=sig+sig-1;

  return(4*Pi^1.25*exp((1-4*MIN_T0*MIN_T0)/(8*h^2)-Pi*A/2)*(1+1/(Pi*A))+2*zeta(sig)*Pi^(-ts/4)*C(sig,t0,h,0)*exp(ts^2/(8*h^2)-A*Pi*ts/2)*(1+1/(A*Pi*ts)));
}

myerfc(x)=
{
   if(x>1000,return(1.9e-434298),return(erfc(x)));
}

ftwiderr(B,t0,h)=
{
   local(bet,X,Y,Z);

   if(t0<exp(exp(1)),print("t0 too small in f~ error."));
   if(B/2>t0,print("B too large in f~ error."));
   bet=1/6+log(log(t0))/log(t0);
   if(B/2<bet*h*h/t0,print("B too small in f~ error."));
   X=(B/2+t0)^bet*exp(-B^2/2/h^2);
   Y=t0^bet*sqrt(Pi/2)*(myerfc(B/2/h/sqrt(2))-myerfc(t0/h/sqrt(2)));
   Z=2^((bet-1)/2)*h^bet*incgam((bet+1)/2,B^2/8/h^2);
   return(24*(X+2^bet*h/B*(Y+Z)));
}

inter_err1(Ns,A,t0,H)=
{
   local(X,Y,Z,bet);
   if(t0<100,print("I use t0=100 in int_err1."));
   if(t0<=exp(exp(1)),print("t0 too small in inter_err1."));
   if(Ns/A>t0,print("t0 too small in inter_err1."));
   bet=1/6+log(log(t0))/log(t0)-1/4;
   /*print("beta="bet);*/
   if(Ns/A<bet*H*H/t0,print("Ns/A too small in inter_err1."));
   X=A/Ns*(Ns/A+t0)^bet*exp(-Ns*Ns/(2*A^2*H^2));
   /*print("X=",X);*/
/* ignore 2nd incgam term */
   Y=t0^bet*incgam(0,Ns*Ns/(2*A^2*H^2));
   /*print("Y=",Y);*/
   Z=2^(bet/2)*H^bet*incgam(bet/2,MIN_T01*MIN_T01/H^2);
   /*print("Z=",Z);*/
   return(6/Pi*(X+2^(bet-1)*(Y+Z)));
}


/* interpolation error due to truncating sum */
inter_err2(sig,t0,A,H)=
{
   local(ts);

   ts=sig+sig-1;
   return(4*zeta(sig)/ts*Pi^((-3-sig-sig)/4)*C(sig,t0,H,0)*exp(ts^2/(8*H^2)-ts*Pi*A/2)+8*Pi^0.25*exp((1-4*MIN_T01*MIN_T01)/(8*H^2)-Pi*A/2));
}

inter_err(sig,Ns,A,t0,H)=
{
   return(inter_err1(Ns,A,t0,H)+inter_err2(sig,t0,A,H));
}

Fmax(x,sig,t0,h)=
{
   local(two_s);
   two_s=sig+sig-1;
   return(zeta(sig)*Pi^(-two_s/4)*C(sig,t0,h,0)*exp(two_s^2/8/h^2-Pi*x*two_s));
}


printpars()=
{
   print("t0=",t0);
   print("A=",A);
   print("N=",N);
   print("B=",B);
   print("h=",h);
   print("K=",K);
   print("M=",M);
   print("M sigma=",msig);
}

main(A,N,upsam,t0,K,inter_step,Ns,H,M)=
{
   local(A1,h,B,N1);

   /*upsam=32;*/
   /*A=4096/21;*/
   A1=A/upsam;
   if(t0<MIN_T0,print("t0 must be at least ",MIN_T0));
   print("minimum T0 = ",MIN_T0);
   /*N=2^20;*/
   while(N<=2^20,
      N1=N/upsam;
      B=N1/A1;
      print("A out=",A);
      print("A (g)=",A1);
      print("B=",B);
      print("N (g)=",N1);
      print("N out=",N);
      print("t0=",t0);
/*
      h=solve(h=B/100,B,log(ftwiderr(B,t0,h))+200);
*/
      h=176431.0/2048.0;

      print("h=",h);
      print("B/h=",B/h);
      print("f~ error=",ftwiderr(B,t0,h));
      print("f^~ error=",fhattwiderr(99,t0,h,A));
/*
      M=103000;
*/
      print("M=",M);
      print("f^ sum error=",bestfhatsumerr(M,t0,h));
/*
      K=solve(k=1,100,log(tayerr(1/2/B,h,k))+130);
      K=ceil(K);
*/
      print("K=",K);
      print("Taylor error=",tayerr(1/2/B,h,K));
      if(B<=h*sqrt(K),print("B needs to be larger"));
      print("g~ error=",maxgtwiderr(B,h,K));
      print("G~ error=",maxGtwiderr(A1,t0,h,K));
      print("Fmax=",Fmax(N1/2/B,99,t0,h));
      print("******************************************************");
      N=N+N;);

   AI=A/inter_step;
   print("A for Interpolation=",AI);
   /*H= 2089.0/16384.0;*/
   print("H=",H);
   /*Ns=70;*/
   print("Ns=",Ns);
   print("Inter err=",inter_err(3,Ns,AI,t0,H));
   print("******************************************************");
}


main(4096/21,2^20,32,3.062*10^10,44,5,70,2089.0/16384.0,104000);


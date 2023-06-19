/*
   zeta-fft-error.gp

   PARI code to calculate parameters for
   windowed FFT based Zeta calculator
*/

MIN_T01=1000 /* used in inter error calc because MIN_T0 causes underflow */
MIN_T0=1e8
/* x is maximum width of approximation = 1/(2B)
   h
   k
*/

lngerr(m,z)=abs(bernfrac(2*m))*2^m/(2*m*(2*m-1)*abs(z)^(2*m-1))

tayerr(x,h,k,M)=
{
   return((2*sqrt(M)-1)*2^((k+5)/2)*Pi^(k+1/2)*h^(k+1)*x^k/gamma((k+2)/2))
}

/* Don't call this with very large sig */

C1(sig,t0,h,k)=
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

/* as per math comp paper */
CC(sig,t0,h,k)=
{
  local(X,Y);

  Y=2^((6*k+5-sig)/4)*Pi^k*exp((1+2*sqrt(2))/(6*t0))*(sig+1/2)^k*(sig+1/2+t0)^((sig-1)/2)*h;
  X=sum(l=0,(sig-1)/2,binomial((sig-1)/2,l)*2^((k+l-1)/2)*h^(k+l+1)+t0^((sig-1-l-l)/2)*incgam((k+l+1)/2,(sig+1/2)^2/2/h^2));
  return(Y+X*2^((6*k+7-sig)/4)*Pi^((2*k-1)/2)*exp((1+2*sqrt(2))/(6*t0)));
}

/* As per TST 19/5/17 */
C(sig,t0,h,k)=
{
  local(X,Y);

  Y=2^((1-sig)/4)*2^((6*k+5-sig)/4)*Pi^k*exp((1+2*sqrt(2))/(6*t0))*(sig+1/2)^k*(sig+1/2+t0)^((sig-1)/2)*h;
  X=2^((1-sig)/4)*sum(l=0,(sig-1)/2,binomial((sig-1)/2,l)*2^((k+l-1)/2)*h^(k+l+1)+t0^((sig-1-l-l)/2)*incgam((k+l+1)/2,(sig+1/2)^2/2/h^2));
  return(Y+X*2^((6*k+7-sig)/4)*Pi^((2*k-1)/2)*exp(sig^2/t0^2+(1+2*sqrt(2))/(6*t0)));
}


/* returns error/k! */
/* should we use t0-B/2 in place of MIN_T0? */
Gtwiderr(A,sig,t0,h,k)=
{
  local(S,s2,s3);

  s2=sig+sig-1;
  s3=sig+sig+1;

  l=0;
  S=sum(l=0,(sig-1)/2,(1+1/(A*Pi*(2*l+1)))*((4*l+1)^2+t0^2)^(k/2)/l!*exp((4*l+1)^2/(8*h*h)-A*Pi*(4*l+1)/2));
  S=S*2^(k+3)*Pi^(k+1);
  S=myexp(log(S)-MIN_T0*MIN_T0/(2*h^2));

  return((2*(1+1/(A*Pi*s2))*C(sig,t0,h,k)*exp(s3^2/(8*h^2)-A*Pi*s2/2)+S)/gamma(k+1));
}


maxGtwiderr(A,t0,h,K)=
{
   local(mx,tp,sig,ttp);
   mx=0;
   for(k=0,K-1,
      sig=3;
      tp=Gtwiderr(A,3,t0,h,k);
      while((sig < 999),
         sig=sig+20;
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

myexp(x)=if(x<-1000,5.1e-435,exp(x));

/* error taking M terms in sum for f^ */
fhatsumerr(M,sig,t0,h)=
{
   local(s2,resid);

   s2=sig+sig-1;
   resid=2*Pi^1.25*myexp(1/(8*h^2)-MIN_T0*MIN_T0/(2*h^2));

   return(resid+C(sig,t0,h,0)*exp(s2^2/(8*h^2))*Pi^(s2*(-0.25))*M^(1-sig)/(sig-1));
}

bestfhatsumerr(M,t0,h)=
{
   local(sig,best_err,err);
   /*printf("In bestfhatsumerr with M=%d t0=%d h=%d\n",M,t0,h);*/
   sig=3;
   best_err=fhatsumerr(M,sig,t0,h);
   while(sig<2000,
         sig=sig+20;
         err=fhatsumerr(M,sig,t0,h);
         if(err<best_err,best_err=err););
   return(best_err);
}

fhattwiderr(sig,t0,h,A)=
{
  local(ts);
  ts=sig+sig-1;

  return(4*Pi^1.25*myexp((1-4*MIN_T0*MIN_T0)/(8*h^2)-Pi*A/2)*(1+1/(Pi*A))+2*zeta(sig)*Pi^(-ts/4)*C(sig,t0,h,0)*exp(ts^2/(8*h^2)-A*Pi*ts/2)*(1+1/(A*Pi*ts)));
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
   bet=log(log(t0))/log(t0);
   /*print("beta="bet);*/
   if(Ns/A<bet*H*H/t0,print("Ns/A too small in inter_err1."));
   X=A/Ns*(Ns/A+t0)^bet*exp(-Ns*Ns/(2*A^2*H^2));
   /*print("X=",X);*/
/* ignore 2nd incgam term */
   Y=t0^bet*incgam(0,Ns*Ns/(2*A^2*H^2));
   /*print("Y=",Y);*/
   Z=2^(bet/2)*H^bet*incgam(bet/2,MIN_T01*MIN_T01/H^2);
   /*print("Z=",Z);*/
   return(24/Pi*(X+2^(bet-1)*(Y+Z)));
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

lgerr=70;

output_dr(str,val)={local(v,e);if(val<1e-307,printf("#define %s ((double) 1e-307)\n",str),e=floor(log(val)/log(10));v=val*10^(-e);printf("#define %s ((double) %4.2fe%d)\n",str,v+0.01,e))};
output_i(str,val)=printf("#define %s ((int) %d)\n",str,val);
output_d(str,val)={local(v,e);e=floor(log(val)/log(10));v=val*10^(-e);printf("#define %s ((double) %10fe%d)\n",str,v,e)};

main(one_over_A,N,upsam,t0,K,inter_step,Ns,H,M,h)=
{
   local(A,A1,B,N1);
  output_d("T0_MAX",t0);
  output_d("T0_MIN",MIN_T0);
  output_i("UPSAM",upsam);
  output_i("N",N);
  output_d("one_over_A",one_over_A*1.0);
  A=1.0/one_over_A;
  A1=A/upsam;
  if(t0<MIN_T0,print("t0 must be at least ",MIN_T0);return(0));
  N1=N/upsam;
  B=N1/A1;
  output_d("h",h);
  et=1/(2*B);
  output_dr("ftwiderr_d",ftwiderr(B,t0,h));
  output_dr("fhattwiderr_d",fhattwiderr(99,t0,h,A));
  output_i("M",M);
  output_dr("fhatsumerr_d",bestfhatsumerr(M,t0,h));
  output_i("K",K);
  te=tayerr(1/2/B,h,K,M);
  output_dr("tayerr_d",te);
  if(B<=h*sqrt(K),print("B needs to be larger");return(0));
  gt=maxgtwiderr(B,h,K);
  output_dr("gtwiderr_d",gt);
  output_dr("Gtwiderr_d",maxGtwiderr(A1,t0,h,K));
  output_dr("Fmaxerr_d",Fmax(N1/2/B,99,t0,h));
  output_i("INTER_SPACING",inter_step);
  AI=A/inter_step;
  output_d("H",H);
  output_i("Ns",Ns);
  output_dr("intererr_d",inter_err(3,Ns,AI,t0,H));
  printf("//************************************\n");
}



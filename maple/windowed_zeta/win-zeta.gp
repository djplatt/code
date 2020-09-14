b1(n,inter_N,t0,A)=t0+(inter_N+n)/A

b(n,inter_N,t0,A,H)=
{
   local(t);

   t=b1(n,inter_N,t0,A);

   return(t^(1/6)*log(t)*exp(-(inter_N+n)^2/(2*H^2*A^2)));
}

inter_err1(inter_N,t0,A,H)=
{
   local(b0);

   b0=b(0,inter_N,t0,A,H);
   return(12*b0/(1-b(1,inter_N,t0,A,H)/b0));
}

C(sig,t0,h)=2^((7-sig)/4)*sqrt(Pi)*exp(1/(2*sig))*gamma((sigma+1)/2)*sum(m=0,(sig-1)/2,t0^((sigma-1)/2-m)*h^(m+1)*2^((m-1)/2)*gamma((m+1)/2)/gamma(m+1)/gamma((sigma-1)/2-m+1))

inter_err2(t0,A,H,sig)=
{
   local(tsig);

   tsig=sig+sig-1;
   return(zeta(sig)/tsig*Pi^((-3-sig-sig)/4)*C(sig,t0,H)*exp(tsig^2/8/H^2-Pi*A*tsig/2));
tsig=2*Pi^0.25*exp((1-4*t0^2)/8/H^2-Pi*A/2);
}

inter_err(inter_N,t0,A,H,sig)=
{
   local(t1,t2);

   t1=inter_err1(inter_N,t0,A,H);
   print("Interpolation sum error=",t1);
   t2=inter_err2(t0,A,H,sig);
   print("Interpolation integral error=",t2);
   return(t1+t2);
}

a(n,B,t0,h)=
{
   local(t);

   t=(n+n+1)*B/2+t0;
   return(t^(1/6)*log(t)*exp(-(n+n+1)^2*B^2/(8*h^2)));
}

f_twiddle_err(B,t0,h)=
{
   local(a0);
+
   a0=a(0,B,t0,h);
   return(6*a0/(1-a(1,B,t0,h)/a0));
}

sum_err(M,h)=
{
   print("ignoring exp(-t0^2) term");
   return(h*sqrt(2*Pi)*M^(-2*Pi^2*h^2*log(M)+1/2)/(2*Pi^2*h^2*log(M)+1/2));
}

f_hat_twiddle_err(h,sig,A)=
{
   local(tsig);

   tsig=sig+sig-1;
   print("ignoring exp(-t0^2) term");
   return(h*sqrt(2*Pi)*zeta(sig)*exp(tsig^2/(8*h^2)-Pi*A*tsig/2)*(1+1/(Pi*A*tsig)));
}

tay_err(K,h,A,xi_1)=
{
   return(h*sqrt(2*Pi)*exp(Pi^2*xi_1^2/(2*h^2))*(xi_1*A*2)^(-K));
}

t0=3*10^10
print("t0=",t0)
A=4096/21 /*131072/6143*/
H=2089/16384 /*12347/65536*/
inter_sig=5 /*3*/
inter_N=70 /*51*/
print("A=",A)
print("Interpolation H=",H)
print("Interpolation sigma=",inter_sig)
print("Interpolation terms=",inter_N)
print("Total Interpolation error=",inter_err(inter_N,t0,A,H,inter_sig))
print("")

M=10
h=50
print("M=",M)
print("h=",h)
print("sum error=",sum_err(M,h))
print("")

xi_1=18
K=16
print("dummy xi=",xi_1)
print("Taylor Terms=",K)
print("tay_err=",tay_err(K,h,A,xi_1))
print("")

twid_sig=8500
print("f_hat_twiddle sigma=",twid_sig)
print("f_hat_twiddle_err=",f_hat_twiddle_err(h,twid_sig,A))
print("")

N0=2^15
B=N0/A
print("N0=",N0)
print("B=",B)
print("f_twiddle_err=",f_twiddle_err(B,t0,h))
print("")

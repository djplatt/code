h=7/32
B=32/5
N=20
Ierr(q,M,h,B,t0)=2*(q/Pi)^(M/2)*zeta(M+0.5)*exp(M^2/2/h^2-2*Pi*B*M)*h*(t0+h/(2*Pi)^0.5+1+1/2/2^0.5)/M
G(n,h,B,t0,N)=(3/2+t0+(N+n)/2/B)^(9/16)*exp(-(N+n)^2/8/B^2/h^2)/Pi/(N+n)
Eerr(q,M,h,B,t0,N)=G(0,h,B,t0,N)/(1-G(1,h,B,t0,N)/G(0,h,B,t0,N))*Pi^0.5*zeta(9/8)*exp(1/6)*2^1.25*(q/2/Pi)^(5/16)
all_err(q,M,h,B,t0,N)=Ierr(q,M,h,B,t0)+Eerr(q,M,h,B,t0,N)

tay_err(N,del,t)=(N+1)*(t+N-0.5)^N*del^N/(N!*(N-1-del*(t+N+0.5))*(N-0.5));

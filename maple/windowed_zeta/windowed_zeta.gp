X(sig,t0,h,k)=sum(l=0,(sig-1)/2,binomial((sig-1)/2,l)*2^((k+l-1)/2)*h^(k+l+1)*t0^((sig-1)/2-l)*incgam((k+l+1)/2,(sig+0.5)^2/(2*h^2)));

C(sig,t0,h,k)=C0(sig,t0,h,k)+C1(sig,t0,h,k);

Fhatsumerr(J,sig,t0,h)=C(sig,t0,h,0)*exp((2*sig-1)^2/(8*h^2))*Pi^((1-2*sig)/4)*J^(1-sig)/(sig-1);

C0(sig,t0,h,k)=2^((6*k+5-sig)/4)*Pi^k*exp((3+sqrt(2))/t0)*((sig+0.5)^k*(sig+0.5+t0)^((sig-1)/2))*h;

C1(sig,t0,h,k)=2^((6*k+7-sig)/4)*Pi^((2*k-1)/2)*exp((3+sqrt(2))/t0)*X(sig,t0,h,0);
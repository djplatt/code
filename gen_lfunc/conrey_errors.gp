
I4(r,q,H,M,R,t0,B)=16*(r+1)*gamma(1+r/2)*q*H*exp(2*M^2/H^2-2*Pi*M*B+(r+2)*Pi*R/4)*(4*zeta(2*M)*(H+2*M+R+3+t0))^(r+2);

E(q,r,R,t0,B,H,N)={local(B0);B0=zeta(3/2)^(2*r);B0*q*(3*r+1)*2^(3*r+1)*((2*R+3+t0)^(3*r)*(1+B*H+B*t0)+2^(r+2)*(2*r+1)!*B*H^(3*r+1)*(N^2/2/H^2)^(2*r))*exp(-N^2/2/H^2);}

inter_err(q,r,R,t0,B,H,N,M)=E(q,r,R,t0,B,H,N)+I4(r,q,H,M,R,t0,B);



/* See Riesel & Vaughn, Platt & Trudgian Brun's constant

Compute sum |g(n)|n^-s via its Euler product. */

gterm(s,p)={local(s1,s2,s3,p2,p3,tmp);p2=p*p;p3=p2*p;s1=p^(-s);s2=s1*s1;s3=s1*s2;tmp=p3-2.0*p2;return(log(1.0+4.0/(p2-2.0*p)*s1+(3*p+2.0)/tmp*s2+2.0/tmp*s3))}

gp(p,p2,p3)=4.0/(p2-2.0*p);
gp2(p,p2,p3)=(3*p+2.0)/(p3-2.0*p2);
gp3(p,p2,p3)=2.0/(p3-2.0*p2);
H(s,P,err)={s1=-s;s2=-2.0*s;s3=-3.0*s;res=log(1.0+0.75*2^s2+0.25*2^s3);forprime(p=3,P,res+=gterm(s,p));exp(res+err)}
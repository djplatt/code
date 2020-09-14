/* see de-Reyna */
rs_err(t,k)=
{
   local(a);
   a=sqrt(t/2/Pi);
   return(a^(-0.24)*2^0.75/7*gamma((k+1)/2)/(10/11*a)^(k+1));
}

/* Berry's form */
c0(r)=cos(Pi/2*(r^2+3/4))/cos(Pi*r)

sec(x)=1/cos(x);


c1(x)=-1/96/Pi^2*(sec(Pi*x)* (-3* Pi^2* x* cos((Pi* (3/4 + x^2))/2) + Pi^3* x^3* sin((Pi* (3/4 + x^2))/2)) + 3* Pi* sec(Pi* x)* (-(Pi^2* x^2* cos((Pi* (3/4 + x^2))/2)) - Pi* sin((Pi* (3/4 + x^2))/2))* tan(Pi* x) - 3* Pi* x* sin((Pi* (3/4 + x^2))/2)* (Pi^2* sec(Pi* x)^3 + Pi^2* sec(Pi* x)* tan(Pi* x)^2) + cos((Pi* (3/4 + x^2))/2)* (5* Pi^3* sec(Pi* x)^3* tan(Pi* x) + Pi^3* sec(Pi* x)* tan(Pi* x)^3));

/* p^2-p-1/16 form */

C0(p)=cos(2*Pi*(p^2-p-1/16))/cos(2*Pi*p);

C1(p)=-1/96/Pi^2*(sec(2*p*Pi)*(-24 *(-1 + 2*p)*Pi^2 *cos(2*(-1/16 - p + p^2)*Pi) + 8 *(-1 + 2*p)^3 *Pi^3 *sin(2*(-1/16 - p + p^2)*Pi)) + 6 *Pi *sec(2*p*Pi)*(-4 *(-1 + 2*p)^2 *Pi^2 *cos(2*(-1/16 - p + p^2)*Pi) - 4 *Pi *sin(2*(-1/16 - p + p^2)*Pi)) *tan(2*p*Pi) - 6 *(-1 + 2*p)*Pi *sin(2*(-1/16 - p + p^2)*Pi)*(4 *Pi^2 *sec(2*p*Pi)^3 + 4 *Pi^2 *sec(2*p*Pi) *tan(2*p*Pi)^2) + cos(2*(-1/16 - p + p^2)*Pi)*(40 *Pi^3 *sec(2*p*Pi)^3 *tan(2*p*Pi) + 8 *Pi^3 *sec(2*p*Pi) *tan(2*p*Pi)^3));

thet(t)=t/2*(log(t/2/Pi)-1)-Pi/8+1/48/t+7/5760/t^3+31/80640/t^5;
/*
imag(lngamma(0.25+I*t/2))-t/2*log(Pi)
*/
rs(t)=
{
   local(a,N,r,th,res,c_0,c_1);
   a=sqrt(t/2/Pi);
   print("a=",a);
   V=floor(a);
   print("V=",V);
   r=1.0-2*(a-V);
   print("r=",r);
   th=thet(t);
   print("theta=",th);
   res=sum(k=1,V,cos(t*log(k)-th)/sqrt(k))*2;
   print("Result of sum is ",res);
   c_0=(-1)^(V+1)/sqrt(a)*c0(r);
   C_0=C0(a-V)/sqrt(a)*(-1)^(V+1);
   print("c0 term is ",c_0);
   print("C0 term is ",C_0);
   c_1=(-1)^(V+1)/(sqrt(a)^3)*c1(r);
   C_1=C1(a-V)/(sqrt(a)^3)*(-1)^(V+1);
   print("c1 term is ",c_1);
   print("C1 term is ",C_1);
   print("Berry form gives ",res+c_0+c_1);
   print("Other form gives ",res+C_0+C_1);
   return(res+c_0+c_1);
}

my_zeta(t)=rs(t)*exp(thet(t)*I)

f(t)=zeta(0.5+I*t)*exp(I*thet(t));

Z(t)=Z1(t/2/Pi);

Z1(t2pi)=Z2(1-2*frac(t2pi));

Z2(x)=cos(0.5*Pi*(x^2+0.75))/cos(Pi*x)


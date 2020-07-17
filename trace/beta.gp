/*
based on betahat =c*sinc^8(ar)(b^2-r^2
*/

beta0(t)=(1/16)*(a^3*(4832*Pi^2*a^7*b^2-3360*Pi^2*a^5*b^2*t^2+1120*Pi^2*a^3*b^2*t^4-280*Pi^2*a*b^2*t^6+70*Pi^2*b^2*t^7-1680*a^5+3360*a^3*t^2-2100*a*t^4+735*t^5)*(1/(a^8*(302*Pi^2*a^2*b^2-105))));

beta1(t)=(1/16)*(a^3*(4944*Pi^2*a^7*b^2-784*Pi^2*a^6*b^2*t-1008*Pi^2*a^5*b^2*t^2-3920*Pi^2*a^4*b^2*t^3+5040*Pi^2*a^3*b^2*t^4-2352*Pi^2*a^2*b^2*t^5+504*Pi^2*a*b^2*t^6-42*Pi^2*b^2*t^7-504*a^5-5880*a^4*t+15120*a^3*t^2-11760*a^2*t^3+3780*a*t^4-441*t^5)*(1/(a^8*(302*Pi^2*a^2*b^2-105))));

beta2(t)=-(1/16)*(a^3*(2224*Pi^2*a^7*b^2-24304*Pi^2*a^6*b^2*t+38640*Pi^2*a^5*b^2*t^2-27440*Pi^2*a^4*b^2*t^3+10640*Pi^2*a^3*b^2*t^4-2352*Pi^2*a^2*b^2*t^5+280*Pi^2*a*b^2*t^6-14*Pi^2*b^2*t^7+19320*a^5-41160*a^4*t+31920*a^3*t^2-11760*a^2*t^3+2100*a*t^4-147*t^5)*(1/(a^8*(302*Pi^2*a^2*b^2-105))));

beta3(t)=(1/16)*(a^3*(32768*Pi^2*a^7*b^2-57344*Pi^2*a^6*b^2*t+43008*Pi^2*a^5*b^2*t^2-17920*Pi^2*a^4*b^2*t^3+4480*Pi^2*a^3*b^2*t^4-672*Pi^2*a^2*b^2*t^5+56*Pi^2*a*b^2*t^6-2*Pi^2*b^2*t^7+21504*a^5-26880*a^4*t+13440*a^3*t^2-3360*a^2*t^3+420*a*t^4-21*t^5)*(1/(a^8*(302*Pi^2*a^2*b^2-105))));

beta(t)=if(t<a,beta0(t),if(t<2*a,beta1(t),if(t<3*a,beta2(t),if(t<4*a,beta3(t),0))));

hhat(t)=(1-beta(t))/(2*(Pi*t)^2);

L(t)=
{
	local(d,l,Lval,ff,pp);
	d=quaddisc(t^2-4);
	l=round(sqrt((t^2-4)/d));
	Lval=2/sqrt(d)*quadregulator(d)*qfbclassno(d);
	ff=factor(l);
	Lval *= prod(i=1,matsize(ff)[1],
		1+(ff[i,1]-kronecker(d,ff[i,1]))*(ff[i,1]^ff[i,2]-1)/(ff[i,1]-1)
	);
	Lval=Lval/l;
        print("Lval(",d,") = ",Lval);
        Lval;
        
}

H()=
{
   local(res,tterm,t,lv);
   t=3;
   res=0;
   tterm=log((t+sqrt(t^2-4))/2)/Pi;
   while(tterm<a,
      res+=L(t)*beta0(tterm)/(2*(Pi*tterm)^2);
      t++;
      tterm=log((t+sqrt(t^2-4))/2)/Pi;);
   while(tterm<2*a,
      res+=L(t)*beta1(tterm)/(2*(Pi*tterm)^2);
      t++;
      tterm=log((t+sqrt(t^2-4))/2)/Pi;);
   while(tterm<3*a,
      res+=L(t)*beta2(tterm)/(2*(Pi*tterm)^2);
      t++;
      tterm=log((t+sqrt(t^2-4))/2)/Pi;);
   while(tterm<4*a,
      res+=L(t)*beta3(tterm)/(2*(Pi*tterm)^2);
      t++;
      tterm=log((t+sqrt(t^2-4))/2)/Pi;);
   print("Last t used was ",t-1);
   return(res/Pi);
};

H1_pow(i)=
{
   local(res,p,lp,lterm);
   p=2;
   n=p^i;
   lp=log(p);
   lterm=i*lp/Pi;
   while(lterm<a,
      n=p^i;
      res+=lp/n*beta0(lterm)/(2*(Pi*lterm)^2);
      p=nextprime(p+1);
      lp=log(p);
      lterm=i*lp/Pi);
   while(lterm<2*a,
      n=p^i;
      res+=lp/n*beta1(lterm)/(2*(Pi*lterm)^2);
      p=nextprime(p+1);
      lp=log(p);
      lterm=i*lp/Pi);
   while(lterm<3*a,
      n=p^i;
      res+=lp/n*beta2(lterm)/(2*(Pi*lterm)^2);
      p=nextprime(p+1);
      lp=log(p);
      lterm=i*lp/Pi);
   while(lterm<4*a,
      n=p^i;
      res+=lp/n*beta3(lterm)/(2*(Pi*lterm)^2);
      p=nextprime(p+1);
      lp=log(p);
      lterm=i*lp/Pi);
   print("Last p^",i," used was ",precprime(p-1));
   return(res);
}

H1()=
{
   local(res,i,h1);
   i=1;
   res=0;
   h1=H1_pow(1);
   while(h1>0,res+=h1;i++;h1=H1_pow(i));
   return(res/Pi);
}

a=1/16;
b=sqrt(6*Pi^2 -1)/2;
c=630*Pi^2*a^3/(302*Pi^2*a^2*b^2-105);

print("a=",a);
print("b=",b);
h0=H();
print("Sum over class numbers = ",h0);
h1=H1();
print("Sum over prime powers  = ",h1);

print("D=",h0+h1);

/*
P=1/Pi*intnum(t=0,[1],(log(2*Pi)-2*real(psi(1+2*I*t)))*h(t));
print("P=",P);
*/

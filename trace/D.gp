/* uses g(x)=1/2pi int h(r) exp(-I r x) dr */

h(t)=if(t==0,1,(sin(Pi*a*t)/(Pi*a*t))^6);

dh(t)=if(t==0,1,6/t*(sin(Pi*a*t)/(Pi*a*t))^5*(cos(a*Pi*t)-sin(a*Pi*t)));

g(t)=if(t<2*Pi*a,2112*a^5*Pi^5-480*a^3*Pi^3*t^2+60*a*Pi*t^4-10*t^5,if(t<4*Pi*a,1632*a^5*Pi^5+1200*a^4*Pi^4*t-1680*a^3*Pi^3*t^2+600*Pi^2*a^2*t^3-90*a*Pi*t^4+5*t^5,(6*Pi*a-t)^5))/(7680*a^6*Pi^6);

dg(t)=if(t<2*Pi*a,-960*a^3*Pi^3*t+240*a*Pi*t^3-50*t^4,if(t<4*Pi*a,1200*a^4*Pi^4-3360*a^3*Pi^3*t+1800*Pi^2*a^2*t^2-360*a*Pi*t^3+25*t^4,-5*(6*Pi*a-t)^4))/(7680*a^6*Pi^6);


L(t)=
{
	local(d,l,Lval,ff);
	d=quaddisc(t^2-4);
	l=round(sqrt((t^2-4)/d));
	Lval=2/sqrt(d)*quadregulator(d)*qfbclassno(d);
        /*print("L(1,chi) for ",d," was ",Lval);*/
	ff=factor(l);
	Lval *= prod(i=1,matsize(ff)[1],
		1+(ff[i,1]-kronecker(d,ff[i,1]))*(ff[i,1]^ff[i,2]-1)/(ff[i,1]-1)
	);
	Lval/l
}

H(maxt)=
{
   local(res,tterm,t,lv);
   t=3;
   res=0;
   tterm=2*log((t+sqrt(t^2-4))/2);
   while(tterm<maxt,
      res+=L(t)*g(tterm);
      t++;
      tterm=2*log((t+sqrt(t^2-4))/2));
   print("Last t used in class numbers was ",t-1);
   return(2*res);
};

H1_pow(i,maxt)=
{
   local(res,p,lp,lterm);
   res=0;
   p=2;
   n=p^i;
   lp=log(p);
   lterm=2*i*lp;
   while(lterm<maxt,
      res+=lp/n*g(lterm);
      p=nextprime(p+1);
      lp=log(p);
      lterm=2*i*lp);
   print("Last p^",i," used in prime sum was ",precprime(p-1));
   return(res);
}

H1(maxt)=
{
   local(res,i,h1);
   i=1;
   res=0;
   h1=H1_pow(1,maxt);
   while(h1>0,res+=h1;i++;h1=H1_pow(i,maxt));
   return(res);
}

a=31/32.;

print("a=",a);
ha=H(6*Pi*a);
print("Sum over class numbers = ",ha);
hb=H1(6*Pi*a);
print("Sum over prime powers  = ",hb);
print("D=",ha+hb);

/*
CC=intnum(t=0,6*Pi*a,hhat(t)*(cosh(t/2)-1));
print("C=",CC);
*/

CC=h(I/2)-h(0);
print("C=",CC);


II=-intnum(t=1e-20,6*Pi*a,dg(t)/sinh(t/2))/6;
print("I=",II);
print("or ",1/6*sum(n=0,10,intnum(r=10^n-1,10^(n+1)-1,r*h(r)*tanh(Pi*r))));

EE=2*intnum(t=0,6*Pi*a,(1/(8*cosh(t/2))+2*cosh(t/2)/(3+6*cosh(t)))*g(t));
print("E= ",EE);
print("or ",2*intnum(r=0,10^10,(1/8+1/(3*sqrt(3))*cosh(Pi*r/3))/cosh(Pi*r)*h(r)));

PP=1/Pi*sum(n=0,10,intnum(t=10^n-1,10^(n+1)-1,(log(2*Pi)-2*real(psi(1+2*I*t)))*h(t)));
print("P= ",PP);
print("or ",g(0)*(log(Pi/2)+2*Euler)-2*intnum(t=1e-10,6*Pi*a,log(4*sinh(t/4))*dg(t))-h(0)/4);


print("I+E+P-h(0)+H+H1-C=",II+EE+PP-h(0)-CC+ha+hb);

read("evs.gp");

tr=sum(n=1,length(evlist_left),h(evlist_left[n]));

print("trace for r<178 = ",tr);

N(t)=t^2/12-2*t/Pi*log(t/exp(1)/sqrt(Pi/2))-131/144;

t0=178;print("estimate for r<infty = ",tr+sum(n=0,1000,t0=t0*1.01;t1=t0*1.01;(N(t1)-N(t0))*h(t0)));print("t0=",t0);
